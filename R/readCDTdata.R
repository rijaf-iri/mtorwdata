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
        stop("Station data is not in a standard unambiguous CDT format")
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
