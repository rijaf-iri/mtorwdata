
ncInfo.dates.table <- function(ncdf, dates.table){
    dates <- dates.table[, -ncol(dates.table), drop = FALSE]
    ldates <- split(dates, col(dates))
    dates <- do.call(paste0, ldates)

    ncfiles <- do.call(sprintf, c(list(fmt = ncdf$format), ldates))

    nc.path <- file.path(ncdf$dir, ncfiles)
    nc.exist <- file.exists(nc.path)
    if(!any(nc.exist)) return(NULL)

    list(dates = dates, ncfiles = nc.path, exist = nc.exist)
}

ncInfo.with.date.range <- function(ncdf, date.range, tstep, minhour = NA){
    dates <- table.format.date.time(tstep, date.range, minhour)
    ncInfo.dates.table(ncdf, dates)
}
