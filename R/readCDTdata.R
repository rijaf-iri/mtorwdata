#' Read CDT station format file.
#'
#' Read CDT station format file.
#' 
#' @param file Full path to the file containing the CDT data
#' @param sep The column's separator of the data
#' @param missing The missing values flag
#'  
#' @return A list object
#' \itemize{
#'   \item{\strong{id}: }{Vector of the points/stations id}
#'   \item{\strong{lon}: }{Vector of the points/stations longitude}
#'   \item{\strong{lat}: }{Vector of the points/stations latitude}
#'   \item{\strong{dates}: }{Vector of the dates or times of the data}
#'   \item{\strong{data}: }{Matrix of the data, row indicates the dates and column the stations}
#' }
#' 
#' @export

readCDTStationFormat <- function(file, sep = ",", missing = "-99"){
    donne <- read.table(file, sep = sep, na.strings = missing,
                        colClasses = "character", stringsAsFactors = FALSE)

    ###############
    seph <- rle(!grepl('[^[:digit:]]', as.character(donne[, 1])))
    ipos <- which(!seph$values & seph$lengths >= 3 & seph$lengths <= 4)
    if(length(ipos) == 0 | ipos[1] != 1){
        stop("Station data is not in a standard CDT format")
    }
    pos <- seph$lengths[ipos[1]]

    ###############
    dat <- list(id = as.character(donne[1, -1]),
                lon = as.numeric(donne[2, -1]),
                lat = as.numeric(donne[3, -1]),
                elv = if(pos == 4) as.numeric(donne[4, -1]) else NULL,
                dates = as.character(donne[-(1:pos), 1]),
                data = local({
                            tmp <- donne[-(1:pos), -1, drop = FALSE]
                            ntmp <- dim(tmp)
                            tmp <- as.numeric(unlist(tmp))
                            dim(tmp) <- ntmp
                            tmp
                    })
                )
    dimnames(dat$data)[[2]] <- NULL
    dat$data <- convert_data_type(dat$data, as.numeric)

    return(dat)
}

#' Filter CDT stations data.
#'
#' Filter two CDT stations data to match the stations and dates.
#' 
#' @param stn_data,grd_data A named list containing the CDT stations data, output from \code{readCDTStationFormat}
#' 
#' @return A list object
#' \itemize{
#'   \item{\strong{stn}: }{The filtered \code{stn_data}}
#'   \item{\strong{grd}: }{The filtered \code{grd_data}}
#' }
#' 
#' @export

filterCDTStationsData <- function(stn_data, grd_data){
    it <- match(grd_data$dates, stn_data$dates)
    it <- it[!is.na(it)]
    if(length(it) == 0)
        stop("Dates do not overlap")

    stn_data$dates <- stn_data$dates[it]
    stn_data$data <- stn_data$data[it, , drop = FALSE]

    it1 <- grd_data$dates %in% stn_data$dates
    grd_data$dates <- grd_data$dates[it1]
    grd_data$data <- grd_data$data[it1, , drop = FALSE]

    id <- match(grd_data$id, stn_data$id)
    id <- id[!is.na(id)]
    if(length(id) == 0)
        stop("Stations do not overlap")

    stn_data$id <- stn_data$id[id]
    stn_data$lon <- stn_data$lon[id]
    stn_data$lat <- stn_data$lat[id]
    stn_data$data <- stn_data$data[, id, drop = FALSE]

    id1 <- grd_data$id %in% stn_data$id
    grd_data$id <- grd_data$id[id1]
    grd_data$lon <- grd_data$lon[id1]
    grd_data$lat <- grd_data$lat[id1]
    grd_data$data <- grd_data$data[, id1, drop = FALSE]

    list(stn = stn_data, grd = grd_data)
}

