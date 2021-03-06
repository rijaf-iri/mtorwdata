#' Precipitation aggregation to hourly
#'
#' Aggregate precipitation accumulation from single scan to hourly.
#' 
#' @param start_time,end_time start and end time of the period to aggregate
#'                   same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param dir5MIN full path to the directory containing the single scan NetCDF files
#' @param dirOUT full path to the directory to save the aggregate data
#' @param min_frac minimum fraction of non missing values
#' @param time_zone the time zone of \code{start_time}, \code{end_time},
#'                  the input and output netCDF files.
#'                  Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' @param in_file_prefix prefix of the input file names. Default "qpe_"
#' @param out_file_prefix prefix of the output file names. Default "precip_"
#' 
#' @export

precip_5min2hour <- function(start_time, end_time,
                             dir5MIN, dirOUT,
                             min_frac = 0.75,
                             time_zone = "Africa/Kigali",
                             in_file_prefix = "qpe_",
                             out_file_prefix = "precip_")
{
    start <- strptime(start_time, "%Y-%m-%d %H:%M", tz = time_zone)
    end <- strptime(end_time, "%Y-%m-%d %H:%M", tz = time_zone) + 43200

    daty <- format(seq(start, end, "day"), "%Y%m%d")

    listFiles <- lapply(daty, function(tt){
        pattern <- paste0(in_file_prefix, tt, ".+\\.nc$")
        list.files(dir5MIN, pattern)
    })
    listFiles <- do.call(c, listFiles)
    if(length(listFiles) == 0) return(NULL)

    # daty <- substr(listFiles, 5, 18)
    daty <- gsub(in_file_prefix, "", listFiles)
    daty <- substr(daty, 1, 14)
    index <- split(seq_along(daty), substr(daty, 1, 10))
    len <- sapply(index, length)

    ifull <- len/12 >= min_frac
    if(all(!ifull)) return(NULL)

    odaty <- names(index)
    ohours <- as.numeric(strptime(odaty, "%Y%m%d%H", tz = time_zone))/3600
    timeUnit <- "hours since 1970-01-01 00:00:00"

    nc <- ncdf4::nc_open(file.path(dir5MIN, listFiles[1]))
    x <- nc$dim$lon$val
    y <- nc$dim$lat$val
    ncdf4::nc_close(nc)
    nx <- length(x)
    ny <- length(y)

    lon <- ncdf4::ncdim_def("lon", "degrees_east", x, longname = "Longitude")
    lat <- ncdf4::ncdim_def("lat", "degrees_north", y, longname = "Latitude")

    for(ix in seq_along(index)){
        if(!ifull[ix]) next
        precip <- lapply(index[[ix]], function(j){
            nc <- ncdf4::nc_open(file.path(dir5MIN, listFiles[j]))
            z <- ncdf4::ncvar_get(nc, 'precip')
            ncdf4::nc_close(nc)
            z[is.na(z)] <- 0
            c(z)
        })

        precip <- do.call(rbind, precip)
        precip <- colSums(precip)
        precip[is.na(precip)] <- -999
        dim(precip) <- c(nx, ny, 1)

        time <- ncdf4::ncdim_def("time", timeUnit, ohours[ix], unlim = TRUE,
                                 calendar = "standard", longname = "Time")
        grd.ncout <- ncdf4::ncvar_def('precip', 'mm', list(lon, lat, time), -999, prec = 'float',
                                      longname = 'Radar estimated precipitation accumulation',
                                      compression = 6)
        out.ncfiles <- file.path(dirOUT, paste0(out_file_prefix, odaty[ix], ".nc"))
        ncout <- ncdf4::nc_create(out.ncfiles, grd.ncout)
        ncdf4::ncvar_put(ncout, grd.ncout, precip)

        ncdf4::ncatt_put(ncout, "lon", "standard_name", "longitude")
        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "standard_name", "latitude")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "description", 'Quantitative precipitation estimation')
        ncdf4::nc_close(ncout)
    }

    return(0)
}

#' Precipitation aggregation from hourly to daily
#'
#' Aggregate hourly precipitation to daily.
#' 
#' @param start_hour,end_hour start and end time of the period to aggregate 
#'                   same time zone as \code{time_zone}, format "YYYY-mm-dd HH:00"
#' @param dirHour full path to the directory containing the NetCDF files
#' @param inFormat format of the input file name. Ex: "precip_adj_\%s\%s\%s\%s.nc"
#' @param dirOUT full path to the directory to save the aggregate data
#' @param outFormat format of the output file name. Ex: "precip_adj_\%s\%s\%s.nc"
#' @param obs_hour observation hour in local time. Default 8h for Rwanda
#' @param min_frac minimum fraction of non missing values
#' @param time_zone the time zone of \code{start_time}, \code{end_time} and the input hourly netCDF files.
#'                 Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' @param long_name long name for the NetCDF variable precip
#' 
#' @export

precip_hour2daily <- function(start_hour, end_hour,
                              dirHour, inFormat,
                              dirOUT, outFormat,
                              obs_hour = 8, min_frac = 0.9,
                              time_zone = "Africa/Kigali",
                              long_name = 'Merged AWS-Radar daily rainfall'
                            )
{
    formatT <- "%Y-%m-%d %H:%M"
    frmt <- c('%Y', '%m', '%d', '%H')
    label <- c('year', 'mon', 'day', 'hour')

    start <- strptime(start_hour, formatT, tz = time_zone)
    end <- strptime(end_hour, formatT, tz = time_zone)

    if(time_zone == "UTC"){
        obs_hour <- obs_hour - 2
    }

    start <- lapply(frmt, function(f) format(start, f))
    end <- lapply(frmt, function(f) format(end, f))
    names(start) <- paste('start', label, sep = '.')
    names(end) <- paste('end', label, sep = '.')
    data.range <- c(start, end)

    ncdf <- list(dir = dirHour, format = inFormat)

    ncinfo <- ncInfo.with.date.range(ncdf, data.range, "hourly", 1)
    if(!any(ncinfo$exist)) return(NULL)

    dates <- ncinfo$dates[ncinfo$exist]
    ncPATH <- ncinfo$ncfiles[ncinfo$exist]

    index <- index.minhr2daily(dates, "hourly", obs_hour)
    len <- sapply(index, length)

    ifull <- len/24 >= min_frac
    if(all(!ifull)) return(NULL)

    odaty <- names(index)
    odays <- as.numeric(as.Date(odaty, "%Y%m%d"))
    timeUnit <- "days since 1970-01-01"

    nc <- ncdf4::nc_open(ncPATH[1])
    x <- nc$dim$lon$val
    y <- nc$dim$lat$val
    ncdf4::nc_close(nc)
    nx <- length(x)
    ny <- length(y)

    lon <- ncdf4::ncdim_def("lon", "degrees_east", x, longname = "Longitude")
    lat <- ncdf4::ncdim_def("lat", "degrees_north", y, longname = "Latitude")

    for(ix in seq_along(index)){
        if(!ifull[ix]) next

        precip <- lapply(index[[ix]], function(j){
            nc <- ncdf4::nc_open(ncPATH[j])
            z <- ncdf4::ncvar_get(nc, 'precip')
            ncdf4::nc_close(nc)
            c(z)
        })
        precip <- do.call(rbind, precip)

        miss <- colSums(!is.na(precip))/24 < min_frac
        if(all(miss)) next

        precip <- colSums(precip)
        precip[miss] <- -999
        precip[is.na(precip)] <- -999
        dim(precip) <- c(nx, ny, 1)

        time <- ncdf4::ncdim_def("time", timeUnit, odays[ix], unlim = TRUE,
                                 calendar = "standard", longname = "Time")
        grd.ncout <- ncdf4::ncvar_def('precip', 'mm', list(lon, lat, time), -999, prec = 'float',
                                      longname = long_name, compression = 6)
        outfrmt <- sprintf(outFormat,
                           substr(odaty[ix], 1, 4),
                           substr(odaty[ix], 5, 6),
                           substr(odaty[ix], 7, 8))
        out.ncfiles <- file.path(dirOUT, outfrmt)
        ncout <- ncdf4::nc_create(out.ncfiles, grd.ncout)
        ncdf4::ncvar_put(ncout, grd.ncout, precip)

        ncdf4::ncatt_put(ncout, "lon", "standard_name", "longitude")
        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "standard_name", "latitude")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "description", long_name)
        ncdf4::nc_close(ncout)
    }

    return(0)
}

#' Precipitation aggregation from daily to monthly
#'
#' Aggregate daily precipitation to monthly.
#' 
#' @param start_day,end_day start and end date of the period to aggregate with format YYYY-mm-dd
#' @param dirDay full path to the directory containing the NetCDF files
#' @param inFormat format of the input file name. Ex: "precip_mrg_\%s\%s\%s.nc"
#' @param dirOUT full path to the directory to save the aggregate data
#' @param outFormat format of the output file name. Ex: "precip_mrg_\%s\%s.nc"
#' @param min_frac minimum fraction of non missing values
#' @param long_name long name for the NetCDF variable precip
#' 
#' @export

precip_daily2month <- function(start_day, end_day,
                               dirDay, inFormat,
                               dirOUT, outFormat,
                               min_frac = 0.9,
                               long_name = 'Merged AWS-Radar monthly rainfall')
{
    frmt <- c('%Y', '%m', '%d')
    label <- c('year', 'mon', 'day')

    start <- as.Date(start_day, "%Y-%m-%d")
    end <- as.Date(end_day, "%Y-%m-%d")
    start <- lapply(frmt, function(f) format(start, f))
    end <- lapply(frmt, function(f) format(end, f))
    names(start) <- paste('start', label, sep = '.')
    names(end) <- paste('end', label, sep = '.')
    data.range <- c(start, end)

    ncdf <- list(dir = dirDay, format = inFormat)

    ncinfo <- ncInfo.with.date.range(ncdf, data.range, "daily")
    if(!any(ncinfo$exist)) return(NULL)

    dates <- ncinfo$dates[ncinfo$exist]
    ncPATH <- ncinfo$ncfiles[ncinfo$exist]

    index <- split(seq_along(dates), substr(dates, 1, 6))
    len <- sapply(index, length)

    odaty <- names(index)
    nbday <- nb.Day.Of.Month(odaty)

    ifull <- len/nbday >= min_frac
    if(all(!ifull)) return(NULL)

    odays <- as.numeric(as.Date(paste0(odaty, "15"), "%Y%m%d"))
    timeUnit <- "days since 1970-01-01"

    nc <- ncdf4::nc_open(ncPATH[1])
    x <- nc$dim$lon$val
    y <- nc$dim$lat$val
    ncdf4::nc_close(nc)
    nx <- length(x)
    ny <- length(y)

    lon <- ncdf4::ncdim_def("lon", "degrees_east", x, longname = "Longitude")
    lat <- ncdf4::ncdim_def("lat", "degrees_north", y, longname = "Latitude")

    for(ix in seq_along(index)){
        if(!ifull[ix]) next

        precip <- lapply(index[[ix]], function(j){
            nc <- ncdf4::nc_open(ncPATH[j])
            z <- ncdf4::ncvar_get(nc, 'precip')
            ncdf4::nc_close(nc)
            c(z)
        })
        precip <- do.call(rbind, precip)

        miss <- colSums(!is.na(precip))/nbday[ix] < min_frac
        if(all(miss)) next

        precip <- colSums(precip)
        precip[miss] <- -999
        precip[is.na(precip)] <- -999
        dim(precip) <- c(nx, ny, 1)

        time <- ncdf4::ncdim_def("time", timeUnit, odays[ix], unlim = TRUE,
                                 calendar = "standard", longname = "Time")
        grd.ncout <- ncdf4::ncvar_def('precip', 'mm', list(lon, lat, time), -999, prec = 'float',
                                      longname = long_name, compression = 6)
        outfrmt <- sprintf(outFormat,
                           substr(odaty[ix], 1, 4),
                           substr(odaty[ix], 5, 6))
        out.ncfiles <- file.path(dirOUT, outfrmt)
        ncout <- ncdf4::nc_create(out.ncfiles, grd.ncout)
        ncdf4::ncvar_put(ncout, grd.ncout, precip)

        ncdf4::ncatt_put(ncout, "lon", "standard_name", "longitude")
        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "standard_name", "latitude")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "description", long_name)
        ncdf4::nc_close(ncout)
    }

    return(0)
}


