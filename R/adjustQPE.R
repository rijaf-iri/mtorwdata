#' Adjust hourly QPE.
#'
#' Adjust hourly QPE using AWS data.
#' 
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param aws_data aws data obtained from \code{awsGetHourlyPrecip}
#'                 A named list of coordinates (list name "coords"), dates (list name "date") and the data (list name "data")
#'                 coordinates: a data.frame with column names "id", "longitude" and "latitude"
#'                 dates: a vector of dates in the format "YYYYmmddHH" in local time
#'                 data: a matrix with row number equals to the length of dates and column number equals to the length of coordinates
#'                 Ex: list(coords = data.frame(id , longitude, latitude), date = vector, data = matrix(nrow = length(date), ncol = nrow(coords)))
#' @param qpe_data named list, directory containing the input netCDF files and the format of the netCDT file names
#'                 Ex: list(dir = "directory/full/path", format = "precip_\%s\%s\%s\%s.nc")
#' @param qpe_adjust named list, directory to save the adjusted QPE and the format of the netCDT file names
#'                   Ex: list(dir = "directory/full/path", format = "precip_adj_\%s\%s\%s\%s.nc")
#' @param pars_adjust named list of the method to be used (list name "method") and other parameters for the adjustment and interpolation (list name "pars").
#'                    The available methods are "Additive", "Multiplicative", "Mixed", "MeanFieldBias" or "KED".
#'                    "Additive": adjustment using an additive error model.
#'                                Default list(method = "Additive").
#'                                \code{pars} can be omitted or list of arguments
#'                                to be passed to the function \code{krige} of the package \code{gstat}
#'                                Ex: list(method = "Additive", pars = list(nmin = 3, nmax = 8))
#'                    "Multiplicative": adjustment using a multiplicative error model.
#'                                \code{pars} can be omitted or list of arguments
#'                                to be passed to the function \code{krige} of the package \code{gstat}
#'                                 Ex: list(method = "Multiplicative", pars = list(nmin = 3, nmax = 8))
#'                    "Mixed": adjustment using a mixed (additive and multiplicative) error model.
#'                                \code{pars} can be omitted or list of arguments
#'                                to be passed to the function \code{krige} of the package \code{gstat}
#'                                 Ex: list(method = "Mixed", pars = list(nmin = 3, nmax = 8))
#'                    "MeanFieldBias": adjustment using one correction factor for the entire domain (Mean Field Bias correction).
#'                                     If \code{pars} is omitted the default is list(method = "linear", minslope = 0.1, minr = 0.5, maxp = 0.01)
#'                                     \code{pars} has the following items:
#'                                     "method": the method used to compute the mean bias fields. Options are "linear", "mean", "median"
#'                                     When using method = "linear" three additional parameters are needed.
#'                                     "minslope": minimum allowable slope
#'                                     "minr": minimum allowable correlation
#'                                     "maxp": maximim allowable p-value
#'                                  Ex: list(method = "MeanFieldBias", pars = list(method = "median"))
#'                  "KED": adjustment using a Kriging with external drift.
#'                       \code{pars} can be omitted or contains a vector candidates of variogram to be fitted to the data, list(models = c("Sph", "Exp", "Gau"))
#'                      Ex: list(method = "KED", pars = list(models = c("Sph", "Exp", "Gau"), nmin = 3, nmax = 8))
#' @param padxy  vector of length 2 representing the number of pixels to be extracted, then aggregated, to get the value of the target pixel.
#'              first element: number of pixels to the left and to the right
#'              second element: number of pixels above and and below
#' @param fun_sp_aggr character, function to be used to aggregate the values of the matched pixels from \code{padxy}.
#'                      Ex: "mean" or "median". Default "median".
#' @param min_aws minimum number of AWS. If the number of AWS with non missing values is less than \code{min_aws},
#'                 no adjustment will be performed.
#' @param min_val minimum value. Only the value greater or equal to \code{min_val} will be used.
#'               For the \code{Additive} and \code{KED} can be 0, otherwise it must be greater than 0
#' @param time_zone the time zone of \code{start_time}, \code{end_time}, the input and output QPE.
#'                  Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' @param ncInfo named list, order of the longitude and latitude dimension in the input netCDF data and the name of the variable.
#'              Default list(ilon = 1, ilat = 2, varid = "precip")
#' 
#' @export

adjustQPE <- function(start_time, end_time, aws_data, qpe_data, qpe_adjust,
                      pars_adjust = list(method = "Additive"),
                      padxy = c(2, 2), fun_sp_aggr = "median",
                      min_aws = 5, min_val = 0.1,
                      time_zone = "Africa/Kigali",
                      ncInfo = list(ilon = 1, ilat = 2, varid = "precip")
                     )
{
    formatT <- "%Y-%m-%d %H:%M"
    formatS <- "%Y%m%d%H"

    start <- strptime(start_time, formatT, tz = time_zone)
    end <- strptime(end_time, formatT, tz = time_zone)
    seqTime <- seq(start, end, "hour")

    yr <- format(seqTime, "%Y")
    mo <- format(seqTime, "%m")
    dy <- format(seqTime, "%d")
    hr <- format(seqTime, "%H")

    ncfileIn <- sprintf(qpe_data$format, yr, mo, dy, hr)
    ncfileOut <- sprintf(qpe_adjust$format, yr, mo, dy, hr)
    ncpathIn <- file.path(qpe_data$dir, ncfileIn)
    ncpathOut <- file.path(qpe_adjust$dir, ncfileOut)

    ix <- file.exists(ncpathIn)
    if(!any(ix)){
        cat("No input netCDF files found\n")
        return(NULL)
    }

    ncpathIn <- ncpathIn[ix]
    ncpathOut <- ncpathOut[ix]
    seqTime <- seqTime[ix]
    seqDate <- format(seqTime, formatS)

    if(time_zone == "UTC"){
        loc2utc <- char_local2utc_time(aws_data$date, formatS)
        aws_data$date <- format(loc2utc, formatS)
    }

    ## get ncdf coords
    nc <- ncdf4::nc_open(ncpathIn[[1]])
    lon <- nc$var[[ncInfo$varid]]$dim[[ncInfo$ilon]]$vals
    lat <- nc$var[[ncInfo$varid]]$dim[[ncInfo$ilat]]$vals
    ncdf4::nc_close(nc)

    nlon <- length(lon)
    nlat <- length(lat)

    ncGRD <- defSpatialPixels(list(lon = lon, lat = lat))
    ijGRD <- getNeighboursIndex(ncGRD, aws_data$coords, padxy)

    if(is.null(ijGRD)){
        cat("AWS and QPE grid do not overlap\n")
        return(NULL)
    }

    ###################

    nclon <- ncdf4::ncdim_def("lon", "degrees_east", lon, longname = "Longitude")
    nclat <- ncdf4::ncdim_def("lat", "degrees_north", lat, longname = "Latitude")

    ###################

    for(jj in seq_along(ncpathIn)){
        ## get AWS data
        it <- which(aws_data$date == seqDate[jj])
        if(length(it) == 0) next
        aws <- aws_data$data[it, ]
        if(all(is.na(aws))) next

        ## get netcdf
        nc <- ncdf4::nc_open(ncpathIn[[jj]])
        rr <- ncdf4::ncvar_get(nc, ncInfo$varid)
        ncdf4::nc_close(nc)

        ## extract netcdf
        qpeSpts <- extractQPEatAWS(ijGRD, rr, aws, fun_sp_aggr)
        ix <- !is.na(qpeSpts$aws) & !is.na(qpeSpts$qpe)
        qpeSpts <- qpeSpts[ix, ]
        im <- qpeSpts$aws >= min_val & qpeSpts$qpe >= min_val
        qpeSpts <- qpeSpts[im, ]
        
        qpe_adj <- rr
        if(length(qpeSpts) >= min_aws){
            qpeGRD <- ncGRD
            qpeGRD$qpe <- c(rr)

            fun_adj <- paste0("qpe_Adjust_", pars_adjust$method)
            fun_adj <- get(fun_adj, mode = "function")
            args <- list(locations = qpeSpts, qpeGRD = qpeGRD)
            if(!is.null(pars_adjust$pars)){
                if(is.list(pars_adjust$pars))
                    args <- c(args, pars_adjust$pars)
            }

            qpe_adj <- do.call(fun_adj, args)
            qpe_adj <- matrix(qpe_adj, nlon, nlat)
        }

        qpe_adj[is.na(qpe_adj)] <- -99
        dim(qpe_adj) <- c(nlon, nlat, 1)

        time <- ncdf4::ncdim_def("time", "seconds since 1970-01-01 00:00:00",
                                 as.numeric(seqTime[jj]), unlim = TRUE,
                                 calendar = "standard", longname = "Time")
        precip <- ncdf4::ncvar_def("precip", "mm", list(nclon, nclat, time), -99,
                                    prec = 'float', compression = 6,
                                    longname = paste("Adjusted QPE using", pars_adjust$method, "method"))

        ncout <- ncdf4::nc_create(ncpathOut[jj], precip)
        ncdf4::ncvar_put(ncout, precip, qpe_adj)

        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "title", "Adjusted Quantitative Precipitation Estimation")
        ncdf4::nc_close(ncout)


        cat(paste("QPE adjustment, time:", seqTime[jj], "done."), "\n")
    }
}


######################

qpe_Adjust_Additive <- function(locations, qpeGRD, ...){
    locations$res <- locations$aws - locations$qpe

    res.grd <- gstat::krige(res~1, locations = locations, newdata = qpeGRD, debug.level = 0, ...)

    out <- qpeGRD$qpe + res.grd$var1.pred
    out[out < 0] <- 0
    return(out)
}

qpe_Adjust_Multiplicative <- function(locations, qpeGRD, ...){
    locations$res <- locations$aws / locations$qpe

    res.grd <- gstat::krige(res~1, locations = locations, newdata = qpeGRD, debug.level = 0, ...)
    res.grd$var1.pred[is.na(res.grd$var1.pred)] <- 1

    out <- qpeGRD$qpe * res.grd$var1.pred
    out[out < 0] <- 0
    return(out)
}

qpe_Adjust_Mixed <- function(locations, qpeGRD, ...){
    locations$eps <- (locations$aws - locations$qpe) / (locations$qpe^2 + 1)
    locations$delta <- ((locations$aws - locations$eps) / locations$qpe) - 1

    eps.grd <- gstat::krige(eps~1, locations = locations, newdata = qpeGRD, debug.level = 0, ...)
    eps.grd$var1.pred[is.na(eps.grd$var1.pred)] <- 0

    delta.grd <- gstat::krige(delta~1, locations = locations, newdata = qpeGRD, debug.level = 0, ...)
    delta.grd$var1.pred[is.na(delta.grd$var1.pred)] <- 0


    out <- (1. + delta.grd$var1.pred) * qpeGRD$qpe + eps.grd$var1.pred
    out[out < 0] <- 0
    return(out)
}

qpe_Adjust_MeanFieldBias <- function(locations, qpeGRD, method = "linear",
                                     minslope = 0.1, minr = 0.5, maxp = 0.01)
{
    ratios <- locations$aws / locations$qpe

    corrfact <- 1
    if(method == "mean") corrfact <- mean(ratios)
    if(method == "median") corrfact <- median(ratios)
    if(method == "linear"){
        linreg <- try(stats::lm(locations$qpe~locations$aws), silent = TRUE)
        if(!inherits(linreg, "try-error")){
            slinr <- summary(linreg)
            pars.lin <- list(slope = slinr$coefficients[2, 1],
                             r = stats::cor(locations$aws, locations$qpe),
                             p = slinr$coefficients[2, 4])
        }else{
            pars.lin <- list(slope = 0, r = 0, p = Inf)
        }

        if(pars.lin$slope > minslope &
           pars.lin$r > minr &
           pars.lin$p < maxp)
        {
            lstsq <- try(stats::lsfit(locations$aws, locations$qpe, intercept = FALSE), silent = TRUE)
            if(!inherits(linreg, "try-error")){
                slope <- lstsq$coefficients
                corrfact <- if(slope == 0) 1 else 1/slope
            }
        }
    }

    out <- corrfact * qpeGRD$qpe
    out[out < 0] <- 0
    return(out)
}

qpe_Adjust_KED <- function(locations, qpeGRD, ...){
    vgm <- NULL
    args <- list(...)
    if("models" %in% names(args)){
        if(!is.null(args$models)){
            if(var(locations$aws - locations$qpe) > 1e-15){
                exp.var <- gstat::variogram(aws ~ qpe, locations = locations, cressie = TRUE)
                vgm <- try(gstat::fit.variogram(exp.var, gstat::vgm(args$models)), silent = TRUE)
                if(inherits(vgm, "try-error")){
                    vgm <- NULL
                }
            }
        }
    }

    default <- c("model", "models", "debug.level", "locations", "newdata")
    args <- args[!names(args) %in% default]
    args <- c(args, list(formula = as.formula("aws ~ qpe"), locations = locations,
                         newdata = qpeGRD, model = vgm, debug.level = 0))
    res.grd <- do.call(gstat::krige, args)

    out <- res.grd$var1.pred
    out[is.infinite(out)] <- NA
    out[out < 0] <- 0
    return(out)
}

