
#' Precipitation accumulation
#'
#' Compute precipitation accumulation for a given duration.
#' 
#' @param tstep time step of the data: "hourly" or "daily" 
#' @param time the time/date at which the total will be computed
#' @param accumul accumulation duration
#' @param min_frac minimum fraction of non missing values
#' @param dirQPE full path to the folder containing the qpe data
#' 
#' @export

qpe_accumulation <- function(tstep, time, accumul, min_frac, dirQPE){
    accumul <- as.numeric(accumul)
    #tstep hourly: hourly, hourly_adj
    #tstep daily: daily, daily_mrg
    if(tstep == "hourly"){
        dirData <- file.path(dirQPE, "hourly")
        temps <- strptime(time, "%Y-%m-%d-%H", tz = "UTC")
        temps <- seq(temps - (accumul - 1) * 3600, temps, 3600)
        temps <- format(temps, "%Y%m%d%H")
        ncfiles <- paste0("precip_", temps, ".nc")
    }else{
        dirData <- file.path(dirQPE, "daily")
        temps <- as.Date(time, "%Y-%m-%d")
        temps <- seq(temps - (accumul - 1), temps, 1)
        temps <- format(temps, "%Y%m%d")
        ncfiles <- paste0("precip_", temps, ".nc")
    }

    ncpath <- file.path(dirData, ncfiles)
    ifile <- file.exists(ncpath)

    if(sum(ifile) / accumul < min_frac) return(NULL)
    ncpath <- ncpath[ifile]

    nc <- ncdf4::nc_open(ncpath[1])
    varid <- nc$var[[1]]$name
    lon <- nc$var[[varid]]$dim[[1]]$vals
    lat <- nc$var[[varid]]$dim[[2]]$vals
    ncdf4::nc_close(nc)

    nx <- length(lon)
    ny <- length(lat)
    miss_data <- matrix(0, nx, ny)
    rr_data <- matrix(0, nx, ny)

    for(jj in seq(accumul)){
        nc <- ncdf4::nc_open(ncpath[jj])
        zval <- ncdf4::ncvar_get(nc, varid)
        ncdf4::nc_close(nc)
        ina <- is.na(zval)
        miss_data <- miss_data + ina
        zval[ina] <- 0
        rr_data <- rr_data + zval
    }

    ina <- miss_data / accumul >= min_frac
    rr_data[ina] <- NA
    out <- list(lon = lon, lat = lat, data = rr_data)

    return(out)
}

#' Precipitation accumulation
#'
#' Write precipitation accumulation to a file in NetCDF format for download.
#' 
#' @param tstep time step of the data: "hourly" or "daily" 
#' @param time the time/date at which the total will be computed
#' @param accumul accumulation duration
#' @param min_frac minimum fraction of non missing values
#' @param dirQPE full path to the folder containing the qpe data
#' @param dirDOWN full path to the folder to write the NetCDF file
#' 
#' @export

down_qpe_accumulation <- function(tstep, time, accumul, min_frac, dirQPE, dirDOWN){
    qpe <- qpe_accumulation(tstep, time, accumul, min_frac, dirQPE)
    if(is.null(qpe)) return(NULL)

    if(tstep == "hourly"){
        temps <- strptime(time, "%Y-%m-%d-%H", tz = "UTC")
        timestamp <- as.numeric(temps)
        timeunit <- "seconds since 1970-01-01 00:00:00"
        temps <- format(temps, "%Y%m%d%H")
        longname <- paste(accumul, "hours Precipitation Accumulation")
    }else{
        temps <- as.Date(time, "%Y-%m-%d")
        timestamp <- as.numeric(temps)
        timeunit <- "days since 1970-01-01"
        temps <- format(temps, "%Y%m%d")
        longname <- paste(accumul, "days Precipitation Accumulation")
    }

    missval <- -99
    ncfile <- paste0("precip_accumul_", temps, ".nc")

    lon <- ncdf4::ncdim_def("lon", "degrees_east", qpe$lon, longname = "Longitude") 
    lat <- ncdf4::ncdim_def("lat", "degrees_north", qpe$lat, longname = "Latitude") 
    time <- ncdf4::ncdim_def("time", timeunit, timestamp, unlim = TRUE,
                             calendar = "standard", longname = "Time")
    grd_ncout <- ncdf4::ncvar_def('precip', 'mm', list(lon, lat, time), missval,
                                  prec = 'float', longname = longname, compression = 9)

    z <- qpe$data
    z[is.na(z)] <- missval
    dim(z) <- c(length(qpe$lon), length(qpe$lat), 1)

    ncpath <- file.path(dirDOWN, ncfile)
    ncout <- ncdf4::nc_create(ncpath, grd_ncout)
    ncdf4::ncvar_put(ncout, grd_ncout, z)

    ncdf4::ncatt_put(ncout, "lon", "axis", "X")
    ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
    ncdf4::ncatt_put(ncout, "time", "axis", "T")
    ncdf4::ncatt_put(ncout, 0, "title", longname)
    ncdf4::nc_close(ncout)

    return(ncfile)
}
