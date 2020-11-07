
Day.Of.Month <- function(year, mon){
    rev((28:31)[!is.na(as.Date(paste(year, mon, 28:31, sep = '-')))])[1]
}

nb.Day.Of.Month <- function(daty){
    nbm <- mapply(Day.Of.Month, substr(daty, 1, 4), substr(daty, 5, 6), USE.NAMES = FALSE)
    as.numeric(nbm)
}

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

index.minhr2daily <- function(times, instep, obs.hour){
    format <- switch(instep, 'minute' = "%Y%m%d%H%M", 'hourly' = "%Y%m%d%H")
    ttn <- as.POSIXct(times, tz = "UTC", format = format)
    ttn <- ttn  - 3600 * obs.hour
    ttn <- format(ttn, format)
    index <- split(seq_along(ttn), substr(ttn, 1, 8))

    return(index)
}

get.range.date.time <- function(date.range, tstep, minhour = NA){
    if(tstep %in% c("daily", "pentad", "dekadal", "monthly")){
        dekpenday <- switch(tstep,
                            "daily" = c("start.day", "end.day"),
                            "pentad" = c("start.pen", "end.pen"),
                            "dekadal" = c("start.dek", "end.dek"),
                            "monthly" = c("start.day", "end.day"))

        start <- date.range[c('start.year', 'start.mon', dekpenday[1])]
        if(tstep == "monthly") start$start.day <- 1
        start <- paste(unlist(start), collapse = "-")
        start <- as.Date(start)

        end <- date.range[c('end.year', 'end.mon', dekpenday[2])]
        if(tstep == "monthly") end$end.day <- 1
        end <- paste(unlist(end), collapse = "-")
        end <- as.Date(end)

        pas <- if(tstep == "monthly") "month" else "day"
    }

    if(tstep == "hourly"){
        start <- date.range[c('start.year', 'start.mon', 'start.day', 'start.hour')]
        if(minhour > 1){
            divh <- start$start.hour %% minhour
            if(divh != 0) start$start.hour <- start$start.hour - divh
        }
        start <- paste(unlist(start), collapse = "-")
        start <- as.POSIXct(start, tz = "UTC", format = "%Y-%m-%d-%H")

        end <- date.range[c('end.year', 'end.mon', 'end.day', 'end.hour')]
        end <- paste(unlist(end), collapse = "-")
        end <- as.POSIXct(end, tz = "UTC", format = "%Y-%m-%d-%H")

        pas <- paste(minhour, "hour")
    }

    if(tstep == "minute"){
        start <- date.range[c('start.year', 'start.mon', 'start.day', 'start.hour', 'start.min')]
        divm <- start$start.min %% minhour
        if(divm != 0) start$start.min <- start$start.min - divm
        start <- paste(unlist(start), collapse = "-")
        start <- as.POSIXct(start, tz = "UTC", format = "%Y-%m-%d-%H-%M")

        end <- date.range[c('end.year', 'end.mon', 'end.day', 'end.hour', 'end.min')]
        end <- paste(unlist(end), collapse = "-")
        end <- as.POSIXct(end, tz = "UTC", format = "%Y-%m-%d-%H-%M")

        pas <- paste(minhour, "min")
    }

    list(start = start, end = end, step = pas)
}

table.format.date.time <- function(tstep, date.range, minhour = NA){
    dates <- get.seq.date.time(date.range, tstep, minhour)

    if(tstep %in% c("daily", "hourly", "minute")){
        doy <- strftime(dates, format = "%j", tz = "UTC")
        format <- switch(tstep,
                         "minute" = "%Y-%m-%d-%H-%M",
                         "hourly" = "%Y-%m-%d-%H",
                         "daily" = "%Y-%m-%d"
                        )
        dates <- format(dates, format)
        dates <- do.call(rbind, strsplit(dates, "-"))
        dates <- cbind(dates, doy)
    }

    if(tstep %in% c("pentad", "dekadal")){
        dates <- format(dates, '%Y-%m-%d')
        dates <- do.call(rbind, strsplit(dates, "-"))

        n <- switch(tstep, "pentad" = 6, "dekadal" = 3)
        xx <- as.numeric(dates[, 3])
        dates <- dates[xx <= n, , drop = FALSE]
        xx <- cbind(stringr::str_pad(rep(1:12, each = n), 2, pad = "0"),
                    stringr::str_pad(rep(1:n, 12), 2, pad = "0"))
        p1 <- paste(dates[, 2], dates[, 3], sep = "-")
        p2 <- paste(xx[, 1], xx[, 2], sep = "-")
        xx <- stringr::str_pad(match(p1, p2), 2, pad = "0")
        dates <- cbind(dates[, 1:2, drop = FALSE], as.numeric(dates[, 3]), xx)
    }

    if(tstep == "monthly"){
        dates <- format(dates, '%Y-%m-%d')
        dates <- do.call(rbind, strsplit(dates, "-"))
    }

    return(dates)
}

get.seq.date.time <- function(date.range, tstep, minhour = NA){
    daty <- get.range.date.time(date.range, tstep, minhour)
    seq(daty$start, daty$end, daty$step)
}

get.range.date.time <- function(date.range, tstep, minhour = NA){
    if(tstep %in% c("daily", "pentad", "dekadal", "monthly")){
        dekpenday <- switch(tstep,
                            "daily" = c("start.day", "end.day"),
                            "pentad" = c("start.pen", "end.pen"),
                            "dekadal" = c("start.dek", "end.dek"),
                            "monthly" = c("start.day", "end.day"))

        start <- date.range[c('start.year', 'start.mon', dekpenday[1])]
        if(tstep == "monthly") start$start.day <- 1
        start <- paste(unlist(start), collapse = "-")
        start <- as.Date(start)

        end <- date.range[c('end.year', 'end.mon', dekpenday[2])]
        if(tstep == "monthly") end$end.day <- 1
        end <- paste(unlist(end), collapse = "-")
        end <- as.Date(end)

        pas <- if(tstep == "monthly") "month" else "day"
    }

    if(tstep == "hourly"){
        start <- date.range[c('start.year', 'start.mon', 'start.day', 'start.hour')]
        if(minhour > 1){
            divh <- start$start.hour %% minhour
            if(divh != 0) start$start.hour <- start$start.hour - divh
        }
        start <- paste(unlist(start), collapse = "-")
        start <- as.POSIXct(start, tz = "UTC", format = "%Y-%m-%d-%H")

        end <- date.range[c('end.year', 'end.mon', 'end.day', 'end.hour')]
        end <- paste(unlist(end), collapse = "-")
        end <- as.POSIXct(end, tz = "UTC", format = "%Y-%m-%d-%H")

        pas <- paste(minhour, "hour")
    }

    if(tstep == "minute"){
        start <- date.range[c('start.year', 'start.mon', 'start.day', 'start.hour', 'start.min')]
        divm <- start$start.min %% minhour
        if(divm != 0) start$start.min <- start$start.min - divm
        start <- paste(unlist(start), collapse = "-")
        start <- as.POSIXct(start, tz = "UTC", format = "%Y-%m-%d-%H-%M")

        end <- date.range[c('end.year', 'end.mon', 'end.day', 'end.hour', 'end.min')]
        end <- paste(unlist(end), collapse = "-")
        end <- as.POSIXct(end, tz = "UTC", format = "%Y-%m-%d-%H-%M")

        pas <- paste(minhour, "min")
    }

    list(start = start, end = end, step = pas)
}
