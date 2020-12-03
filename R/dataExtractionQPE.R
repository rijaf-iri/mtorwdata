#' Extract QPE
#'
#' Extract QPE.
#' 
#' @param dirQPE path to qpe folder
#' @param dirDOWN path to save the extracted qpe
#' @param timestep time step of the data
#' @param timerange time range, list(start, end)
#' @param geom list, geom to extract
#' @param min_frac minimum fraction of non missing values, for pentad and dekadal data
#' 
#' @return a JSON object
#' 
#' @export

extractQPE <- function(dirQPE, dirDOWN, timestep, timerange, geom, min_frac = 1.0)
{
    if(geom$type == 'points'){
        padxy <- do.call(c, geom$padxy)
        padxy <- as.numeric(padxy)
        pts_id <- as.character(do.call(c, geom$id))
        if(geom$extype == "mpoints")
            pts_id <- paste0("Points_", seq_along(pts_id))
        pts_x <- as.numeric(do.call(c, geom$x))
        pts_y <- as.numeric(do.call(c, geom$y))
        pts <- data.frame(id = pts_id, x = pts_x, y = pts_y)

        geomObj <- list(type = geom$type, points = pts, padxy = padxy)
    }else{
        shpf <- rgdal::readOGR(dsn = geom$dsn, layer = geom$layer, verbose = FALSE)

        if(geom$extype == 'ushapefile'){
            if(geom$multishp == 'single'){
                shpf$ExtrAttr <- "VAL"
                attr_values <- "VAL"
                attr_name <- "ExtrAttr"
                nom <- geom$layer
            }
        }else{
            attr_values <- as.character(do.call(c, geom$id))
            attr_name <- geom$field

            iattr <- as.character(shpf@data[, attr_name])
            iname <- as.character(shpf@data[, geom$nameID])
            nom <- iname[match(attr_values, iattr)]
            nom <- trimws(nom)
            nom <- gsub("[^[:alnum:]]", "", nom)
            # nom <- iconv(nom, from = "UTF-8", to = 'ASCII//TRANSLIT')
            nom <- stringi::stri_trans_general(str = nom, id = "Latin-ASCII")
        }

        geomObj <- list(type = geom$type, shp = shpf, attr_values = attr_values,
                        attr_name = attr_name, names = nom, spavg = geom$spavg)
    }

    ncInfo <- switch(timestep,
                     "minute" = c("5minutes", "qpe_%s%s%s%s.+\\.nc$"),
                     "hourly" = c("hourly", "precip_%s%s%s%s.nc"),
                     "daily" = c("daily", "precip_%s%s%s.nc"),
                     "pentad" = c("daily", "precip_%s%s%s.nc"),
                     "dekadal" = c("daily", "precip_%s%s%s.nc"),
                     "monthly" = c("monthly", "precip_%s%s.nc")
                    )
    min_frac <- as.numeric(min_frac)

    ret <- extractQPEnc(timestep, timerange, dirQPE, ncInfo, dirDOWN, geomObj, min_frac)

    if(is.null(ret)){
        out <- "no-data"
        type <- "no-data"
    }else{
        out <- ret
        type <- "csv"
        if(geomObj$type == 'polys'){
            if(geomObj$spavg == 'gridded'){
                out <- "extracted_qpe.zip"
                wd <- getwd()
                setwd(dirDOWN)
                utils::zip(zipfile = out, files = ret)
                setwd(wd)
                type <- "zip"
            }
        }
    }

    convJSON(list(type = type, file = out, dir = dirDOWN))
}

#####################

extractQPEnc <- function(timestep, timerange, dirQPE, ncInfo,
                         dirOUT, geomObj, min_frac = 1.0)
{
    ncdata <- get_qpe_ncdata_files(timestep, timerange, dirQPE, ncInfo)
    if(is.null(ncdata)) return(NULL)

    daty <- ncdata$date
    ncpath <- ncdata$path
    times <- ncdata$time
    nbDay <- ncdata$nbDay

    ncfile <- if(timestep %in% c("pentad", "dekadal")) ncpath[[1]][1] else ncpath[1]
    varInfo <- get_qpe_ncdata_infos(ncfile)

    gridObj <- defSpatialPixels(list(lon = varInfo$lon, lat = varInfo$lat), regrid = TRUE)

    #######

    if(geomObj$type == 'points'){
        ret <- extractGeomPoints(gridObj, geomObj$points, geomObj$padxy)
        if(is.null(ret)) return(NULL)
        headinfo <- ret$headinfo
        ij2xtr <- ret$ij2xtr

        ij2xtr <- lapply(ij2xtr, function(x) x[!is.na(x)])
        nonZero <- sapply(ij2xtr, length) > 0
        ij2xtr <- ij2xtr[nonZero]

        fileout <- paste0("extracted_qpe_", daty[1], "-",
                          daty[length(daty)], ".csv")
        fileCSV <- file.path(dirOUT, fileout)
        out <- cbind(c("Points", "LON", "LAT"), t(headinfo))
        write.table(out, fileCSV, sep = ",", col.names = FALSE,
                                row.names = FALSE, quote = FALSE)
    }else{
        ret <- extractGeomPolys(gridObj, geomObj$shp, geomObj$attr_name, geomObj$attr_values)
        if(is.null(ret)) return(NULL)
        ij2xtr <- ret$ij2xtr

        nonNull <- !sapply(ij2xtr, is.null)
        ij2xtr <- ij2xtr[nonNull]
        ij2xtr <- lapply(ij2xtr, function(x) x[order(x[, "value"]), , drop = FALSE])

        if(geomObj$spavg == "average"){
            fileout <- paste0("extracted_qpe_", daty[1], "-",
                              daty[length(daty)], ".csv")
            fileCSV <- file.path(dirOUT, fileout)
            headinfo <- ret$headinfo
            headinfo[, 1] <- geomObj$names
            out <- cbind(c("Points", "LON", "LAT"), t(headinfo))
            write.table(out, fileCSV, sep = ",", col.names = FALSE,
                                    row.names = FALSE, quote = FALSE)
        }else{
            # ncdf
            headinfo <- lapply(ij2xtr, function(x){
                xx <- gridObj@coords[x[, 'value'], , drop = FALSE]
                nl <- seq(nrow(x))
                reshapeXYZ2Matrix(cbind(xx, nl))
            })

            fileout <- paste0(geomObj$names[nonNull], ".nc")
            ncout <- lapply(seq_along(headinfo), function(ii){
                ncfile <- file.path(dirOUT, fileout[ii])
                lon <- ncdf4::ncdim_def("lon", "degrees_east", headinfo[[ii]]$x,
                                        longname = "Longitude")
                lat <- ncdf4::ncdim_def("lat", "degrees_north", headinfo[[ii]]$y,
                                        longname = "Latitude")

                unit <- if(timestep %in% c("minute", "hourly"))
                            "seconds since 1970-01-01 00:00:00"
                        else
                            "days since 1970-01-01"
                temps <- as.numeric(times)
                time <- ncdf4::ncdim_def("time", unit, temps, unlim = TRUE,
                                         calendar = "standard", longname = "Time")

                grd_nc <- ncdf4::ncvar_def('precip', 'mm', list(lon, lat, time), -99,
                                            longname = varInfo$longname,
                                            prec = 'float', compression = 9)
                nc <- ncdf4::nc_create(ncfile, grd_nc)
                ncdf4::ncatt_put(nc, "lon", "axis", "X")
                ncdf4::ncatt_put(nc, "lat", "axis", "Y")
                ncdf4::ncatt_put(nc, "time", "axis", "T")
                ncdf4::ncatt_put(nc, 0, "title", varInfo$longname)

                list(nc = nc, grd = grd_nc)
            })
        }
    }

    ########

    for(jj in seq_along(ncpath)){
        if(timestep %in% c("pentad", "dekadal")){
            minFrac <- length(ncpath[[jj]])/nbDay[jj] >= min_frac
            if(minFrac){
                ncs <- lapply(ncpath[[jj]], function(pth){
                    nc <- ncdf4::nc_open(pth)
                    val <- ncdf4::ncvar_get(nc, varInfo$name)
                    ncdf4::nc_close(nc)
                    c(val)
                })
                ncs <- do.call(rbind, ncs)
                ina <- colSums(is.na(ncs))
                ina <- ina/nbDay[jj] >= min_frac
                val <- colSums(ncs, na.rm = TRUE)
                val[ina] <- NA
                dim(val) <- c(varInfo$nlon, varInfo$nlat)
            }else{
                val <- matrix(NA, nrow = varInfo$nlon, ncol = varInfo$nlat)
            }
        }else{
            nc <- ncdf4::nc_open(ncpath[jj])
            val <- ncdf4::ncvar_get(nc, varInfo$name)
            ncdf4::nc_close(nc)
        }

        ##########

        if(geomObj$type == 'points'){
            extdat <- sapply(ij2xtr, function(ij){
                mat <- val[ij]
                if(length(ij) > 1)
                    mat <- mean(mat, na.rm = TRUE)
                return(mat)
            })

            if(!all(nonZero)){
                tmp <- rep(NA, length(nonZero))
                tmp[nonZero] <- extdat
                extdat <- tmp
            }
        
            extdat[is.na(extdat)] <- -99
            out <- c(daty[jj], round(extdat, 1))
            out <- paste(out, collapse = ",")
            cat(out, file = fileCSV, sep = "\n", append = TRUE)
        }else{
            extdat <- lapply(seq_along(ij2xtr), function(ii){
                ij <- ij2xtr[[ii]]
                mat <- val[ij[, "value"]]

                if(geomObj$spavg == "average"){
                    if(nrow(ij) > 1){
                        mat <- mat * ij[, "weight"]
                        mat <- mat[!is.na(mat)]
                        mat <- if(length(mat) > 0) sum(mat) else NA
                    }
                }else{
                    mat <- mat[headinfo[[ii]]$z]
                    mat[is.na(mat)] <- -99
                    dim(mat) <- c(dim(headinfo[[ii]]$z), 1)
                }

                return(mat)
            })

            if(geomObj$spavg == "average"){
                extdat <- do.call(c, extdat)

                if(!all(nonNull)){
                    tmp <- rep(NA, length(nonNull))
                    tmp[nonNull] <- extdat
                    extdat <- tmp
                }

                extdat[is.na(extdat)] <- -99
                out <- c(daty[jj], round(extdat, 1))
                out <- paste(out, collapse = ",")
                cat(out, file = fileCSV, sep = "\n", append = TRUE)
            }else{
                for(ii in seq_along(extdat)){
                    ncdf4::ncvar_put(ncout[[ii]]$nc, ncout[[ii]]$grd, extdat[[ii]],
                                     start = c(1, 1, jj), count = c(-1, -1, 1))
                }
            }
        }
    }

    if(geomObj$type == 'polys'){
        if(geomObj$spavg == 'gridded'){
            for(ii in seq_along(ncout))
                ncdf4::nc_close(ncout[[ii]]$nc)
        }
    }

    return(fileout)
}

#####################

get_qpe_ncdata_files <- function(timestep, timerange, dirQPE, ncInfo){
    dirDAT <- file.path(dirQPE, ncInfo[1])

    timeR <- lapply(timerange, strsplit, split = "-")
    timeR <- lapply(timeR, "[[", 1)

    if(timestep == "minute"){
        timeR <- lapply(timeR, function(x){
            t <- paste0(x[1:5], collapse = "-")
            strptime(t, "%Y-%m-%d-%H-%M", tz = "UTC")
        })

        daty <- seq(timeR$start, timeR$end, "hour")
        daty <- lapply(c('%Y', '%m', '%d', '%H'), function(f) format(daty, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]], daty[[3]], daty[[4]])
        ncfiles <- lapply(ncfiles, function(nc){
            list.files(dirDAT, nc)
        })
        ncfiles <- do.call(c, ncfiles)
        if(length(ncfiles) == 0) return(NULL)
        daty <- substr(ncfiles, 5, 18)
        times <- strptime(daty, "%Y%m%d%H%M%S", tz = "UTC")
    }

    if(timestep == "hourly"){
        timeR <- lapply(timeR, function(x){
            t <- paste0(x[1:4], collapse = "-")
            strptime(t, "%Y-%m-%d-%H", tz = "UTC")
        })

        times <- seq(timeR$start, timeR$end, "hour")
        daty <- lapply(c('%Y', '%m', '%d', '%H'), function(f) format(times, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]], daty[[3]], daty[[4]])
    }

    if(timestep == "daily"){
        timeR <- lapply(timeR, function(x){
            t <- paste0(x[1:3], collapse = "-")
            as.Date(t)
        })

        times <- seq(timeR$start, timeR$end, "day")
        daty <- lapply(c('%Y', '%m', '%d'), function(f) format(times, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]], daty[[3]])
    }

    if(timestep == "pentad"){
        start <- c('01', '06', '11', '16', '21', '26')
        timeR$start[3] <- start[as.numeric(timeR$start[3])]
        dmon <- day_of_month(timeR$end[1], timeR$end[2])
        end <- c('05', '10', '15', '20', '25', dmon)
        timeR$end[3] <- end[as.numeric(timeR$end[3])]
        timeR <- lapply(timeR, function(x){
            t <- paste0(x[1:3], collapse = "-")
            as.Date(t)
        })

        times <- seq(timeR$start, timeR$end, "day")
        daty <- lapply(c('%Y', '%m', '%d'), function(f) format(times, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]], daty[[3]])
    }

    if(timestep == "dekadal"){
        start <- c('01', '11', '21')
        timeR$start[3] <- start[as.numeric(timeR$start[3])]
        dmon <- day_of_month(timeR$end[1], timeR$end[2])
        end <- c('10', '20', dmon)
        timeR$end[3] <- end[as.numeric(timeR$end[3])]
        timeR <- lapply(timeR, function(x){
            t <- paste0(x[1:3], collapse = "-")
            as.Date(t)
        })

        times <- seq(timeR$start, timeR$end, "day")
        daty <- lapply(c('%Y', '%m', '%d'), function(f) format(times, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]], daty[[3]])
    }

    if(timestep == "monthly"){
        timeR <- lapply(timeR, function(x){
            t <- paste0(c(x[1:2], 15), collapse = "-")
            as.Date(t)
        })

        times <- seq(timeR$start, timeR$end, "month")
        daty <- lapply(c('%Y', '%m'), function(f) format(times, f))
        ncfiles <- sprintf(ncInfo[2], daty[[1]], daty[[2]])
    }

    ######

    if(timestep == "minute"){
        it <- times >= timeR$start & times <= timeR$end
    }else{
        it <- file.exists(file.path(dirDAT, ncfiles))
    }

    if(!any(it)) return(NULL)
    times <- times[it]
    ncfiles <- ncfiles[it]
    ncpath <- file.path(dirDAT, ncfiles)

    ######

    if(timestep %in% c("pentad", "dekadal")){
        infos <- switch(timestep,
                        "pentad" = list(c(1, 5, 10, 15, 20, 25, 31),
                                        nb_day_of_pentad),
                        "dekadal" = list(c(1, 10, 20, 31),
                                         nb_day_of_dekad)
                       )
        yymm <- format(times, "%Y%m")
        dd <- as.numeric(format(times, '%d'))
        dd <- cut(dd, infos[[1]], labels = FALSE, include.lowest = TRUE)
        index <- split(seq_along(dd), paste0(yymm, dd))
        ncpath <- lapply(index, function(ix) ncpath[ix])
        daty <- names(index)

        nbDay <- infos[[2]](daty)

        ## time for netcdf
        yymm <- substr(daty, 1, 6)
        dd <- as.numeric(substr(daty, 7, 7))
        brks <- infos[[1]]
        brks <- brks[-length(brks)]
        times <- paste0(yymm, '-', brks[dd])
        times <- as.Date(times, format = "%Y%m-%d")
    }else{
        frmt <- switch(timestep,
                        "minute" = '%Y%m%d%H%M',
                        "hourly" = '%Y%m%d%H',
                        "daily" = '%Y%m%d',
                        "monthly" = '%Y%m'
                    )
        ## minute and hourly convert to Kigali time
        if(timestep %in% c("minute", "hourly")){
            daty <- time_utc2local_char(times, frmt)
        }else{
            daty <- format(times, frmt)
        }

        nbDay <- NULL
    }

    list(date = daty, path = ncpath, time = times, nbDay = nbDay)
}

get_qpe_ncdata_infos <- function(ncfile, varid = NULL){
    varInfo <- NULL
    nc <- ncdf4::nc_open(ncfile)
    varInfo$name <- if(is.null(varid)) nc$var[[1]]$name else varid
    varInfo$lon <- nc$var[[varInfo$name]]$dim[[1]]$vals
    varInfo$lat <- nc$var[[varInfo$name]]$dim[[2]]$vals
    varInfo$longname <- nc$var[[varInfo$name]]$longname
    varInfo$units <- nc$var[[varInfo$name]]$units
    varInfo$prec <- nc$var[[varInfo$name]]$prec
    ncdf4::nc_close(nc)

    varInfo$nlon <- length(varInfo$lon)
    varInfo$nlat <- length(varInfo$lat)

    varInfo
}
