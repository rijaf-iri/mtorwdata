
#' Extract CAPPI
#'
#' Extract CAPPI for single scan over a set of points
#' 
#' @param cappi_dir Full path to folder containing the CAPPI netCDF files
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, 
#'                  single scan format "YYYY-mm-dd HH:MM", hourly format "YYYY-mm-dd HH:00"
#' @param fields A vector of fields to extract
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param padxy A vector of the padding to use, in order "lon", "lat".
#'              Default c(0, 0), no padding applied.
#' @param fun_sp Character, function to be used for the padding. Options: "mean", "median", "max", "min"
#'
#'  @return A data.frame
#' 
#' @export

extractCAPPI <- function(cappi_dir, start_time, end_time,
                         fields, points, padxy = c(0, 0),
                         fun_sp = "mean")
{
    start <- strptime(start_time, "%Y-%m-%d %H:%M")
    end <- strptime(end_time, "%Y-%m-%d %H:%M")

    ##########

    timerange <- list(start = format(start, "%Y-%m-%d-%H-%M"),
                      end = format(end, "%Y-%m-%d-%H-%M"))
    timeR <- lapply(timerange, strsplit, split = "-")
    timeR <- lapply(timeR, "[[", 1)
    timeR <- lapply(timeR, function(x){
        t <- paste0(x[1:5], collapse = "-")
        strptime(t, "%Y-%m-%d-%H-%M")
    })


    daty <- seq(timeR$start, timeR$end, "hour")
    daty <- lapply(c('%Y', '%m', '%d', '%H'), function(f) format(daty, f))
    ncfiles <- sprintf("cappi_%s%s%s%s.+\\.nc$", daty[[1]], daty[[2]], daty[[3]], daty[[4]])
    ncfiles <- lapply(ncfiles, function(nc){
        list.files(cappi_dir, nc)
    })

    ncfiles <- do.call(c, ncfiles)
    if(length(ncfiles) == 0){
        cat("No netCDF CAPPI files found\n")
        return(NULL)
    }
    daty <- gsub("[^[:digit:]]", "", ncfiles)
    times <- strptime(daty, "%Y%m%d%H%M%S")

    ##########
    it <- times >= timeR$start & times <= timeR$end

    if(!any(it)) return(NULL)
    times <- times[it]
    ncfiles <- ncfiles[it]
    ncpath <- file.path(cappi_dir, ncfiles)
    ##########

    daty <- format(times, '%Y%m%d%H%M')

    ##########

    nc <- ncdf4::nc_open(ncpath[1])
    var_name <- nc$var[[1]]$name
    lon <- nc$var[[var_name]]$dim[[1]]$vals
    lat <- nc$var[[var_name]]$dim[[2]]$vals
    ncdf4::nc_close(nc)

    gridObj <- defSpatialPixels(list(lon = lon, lat = lat))

    ##########

    fun_sp <- get(fun_sp, mode = "function")

    names(points) <- c('id', 'x', 'y')
    ijGRD <- getNeighboursIndex(gridObj, points, padxy)

    if(is.null(ijGRD)){
        cat("Points and CAPPI grid do not overlap\n")
        return(NULL)
    }

    ij2xtr <- lapply(ijGRD$ij2xtr, function(x) x[!is.na(x)])
    nonZero <- sapply(ij2xtr, length) > 0
    ij2xtr <- ij2xtr[nonZero]

    ##########
    don <- lapply(seq_along(ncpath), function(jj){
        nc <- ncdf4::nc_open(ncpath[jj])
        don <- lapply(fields, function(field) ncdf4::ncvar_get(nc, field))
        ncdf4::nc_close(nc)
        names(don) <- fields

        extdat <- lapply(ij2xtr, function(ij){
            mat <- lapply(don, function(x) x[ij])
            if(length(ij) > 1)
                mat <- lapply(mat, fun_sp, na.rm = TRUE)
            
            do.call(c, mat)
        })
        extdat <- do.call(rbind, extdat)
        extdat[is.nan(extdat) | is.infinite(extdat)] <- NA

        if(!all(nonZero)){
            tmp <- matrix(NA, length(nonZero), length(fields))
            dimnames(tmp)[[2]] <- dimnames(extdat)[[2]]
            tmp[nonZero, ] <- extdat
            extdat <- tmp
        }

        extdat
    })

    don <- lapply(seq_along(don), function(jj){
        cbind(points, time = daty[jj], don[[jj]])
    })

    do.call(rbind, don)
}
