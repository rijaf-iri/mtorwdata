
convJSON <- function(obj, ...){
    args <- list(...)
    if(!'pretty' %in% names(args)) args$pretty <- TRUE
    if(!'auto_unbox' %in% names(args)) args$auto_unbox <- TRUE
    if(!'na' %in% names(args)) args$na <- "null"
    args <- c(list(x = obj), args)
    json <- do.call(jsonlite::toJSON, args)
    return(json)
}

convCSV <- function(obj, col.names = TRUE){
    filename <- tempfile()
    write.table(obj, filename, sep = ",", na = "", col.names = col.names,
                row.names = FALSE, quote = FALSE)
    don <- readLines(filename)
    unlink(filename)
    don <- paste0(don, collapse = "\n")

    return(don)
}

##################

char_utc2local_time <- function(dates, format, tz = "Africa/Kigali"){
    x <- strptime(dates, format, tz = "UTC")
    x <- as.POSIXct(x)
    x <- format(x, format, tz = tz)
    x <- strptime(x, format, tz = tz)
    x
}

time_utc2local_char <- function(dates, format, tz = "Africa/Kigali"){
    x <- as.POSIXct(dates)
    x <- format(x, format, tz = tz)
    x
}

char_local2utc_time <- function(dates, format, tz = "Africa/Kigali"){
    x <- strptime(dates, format, tz = tz)
    x <- as.POSIXct(x)
    x <- format(x, format, tz = "UTC")
    x <- strptime(x, format, tz = "UTC")
    x
}

time_local2utc_char <- function(dates, format){
    x <- as.POSIXct(dates)
    x <- format(x, format, tz = "UTC")
    x
}

##################

day_of_month <- function(year, mon){
    end_mon <- as.Date(paste(year, mon, 28:31, sep = '-'))
    rev((28:31)[!is.na(end_mon)])[1]
}

nb_day_of_month <- function(daty){
    year <- substr(daty, 1, 4)
    mon <- substr(daty, 5, 6)
    nbm <- mapply(day_of_month, year, mon, USE.NAMES = FALSE)
    as.numeric(nbm)
}

nb_day_of_pentad <- function(daty){
    day <- as.numeric(substr(daty, 7, 7))
    nbp <- rep(5, length(daty))
    nbp[day >= 6] <- nb_day_of_month(daty[day == 6]) - 25
    return(nbp)
}

nb_day_of_dekad <- function(daty){
    day <- as.numeric(substr(daty, 7, 7))
    nbd <- rep(10, length(daty))
    nbd[day == 3] <- nb_day_of_month(daty[day == 3]) - 20
    return(nbd)
}

##################

defSpatialPixels <- function(grd_Coords,
                             projCRS = sp::CRS(as.character(NA)),
                             regrid = FALSE)
{
    if(regrid){
        x <- grd_Coords$lon
        xrg <- diff(range(diff(x)))
        if(xrg > 0.0001){
            xr <- range(x)
            x <- seq(xr[1], xr[2], length.out = length(x))
        }
        y <- grd_Coords$lat
        yrg <- diff(range(diff(y)))
        if(yrg > 0.0001){
            yr <- range(y)
            y <- seq(yr[1], yr[2], length.out = length(y))
        }

        grd0 <- expand.grid(lon = x, lat = y)
        sp::coordinates(grd0) <- ~lon+lat
        grd <- sp::SpatialPixels(points = grd0, tolerance = 0.0002, proj4string = projCRS)
    }else{
        grd0 <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
        sp::coordinates(grd0) <- ~lon+lat

        foo <- function(tol)
                    sp::SpatialPixels(points = grd0, tolerance = tol, proj4string = projCRS)
        grd <- try(foo(sqrt(sqrt(.Machine$double.eps))), silent = TRUE)
        if(inherits(grd, "try-error")) grd <- foo(0.005)
    }

    return(grd)
}

##################

reshapeXYZ2Matrix <- function(df){
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    names(df) <- c('x', 'y', 'z')
    x <- sort(unique(df$x))
    y <- sort(unique(df$y))
    z <- reshape2::acast(df, x~y, value.var = "z")
    dimnames(z) <- NULL
    return(list(x = x, y = y, z = z))
}

