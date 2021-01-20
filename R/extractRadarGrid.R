
#' Extract radar Cartesian data.
#'
#' Extract radar Cartesian data over a set of points.
#' 
#' @param url The URL of the server. Ex: "http://192.168.1.10:8080"
#' @param source A partial path of the directory containing the folders dates of the mdv files.
#'              Ex: "radarCart/ops", "ctrec" or "radarCart/echo_tops"
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param fields Vector of the fields to extract.
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param levels A vector of the index of altitudes to be extracted in integer, or -1 to extract all available altitude
#' @param padxyz A vector of the padding to use in number of pixels, in order "lon", "lat", "alt"
#' @param fun_sp Character, function to be used for the padding. Options: "mean", "median", "max", "min"
#' @param time_zone the time zone of \code{start_time}, \code{end_time} and the output QPE.
#'                 Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' 
#' @return A named list
#'        points: the original points used to extract the data
#'        date: list of dates of the extracted data
#'        altitude: list of the altitudes of the extracted data
#'        data: list of the fields in the form of 3d array
#'             dimension: (len(date) x len(altitude) x len(points))
#' 
#' @export

extractRadarGrid <- function(url, source,
                             start_time, end_time,
                             fields, points, levels = -1,
                             padxyz = c(0, 0, 0), fun_sp = "mean",
                             time_zone = "Africa/Kigali")
{
    on.exit(curl::handle_reset(handle))

    #######

    args <- list(source = source,
                 start_time = start_time,
                 end_time = end_time,
                 fields = fields,
                 points = points,
                 levels = levels,
                 padxyz = padxyz,
                 fun_sp = fun_sp,
                 time_zone = time_zone
                )
    args <- jsonlite::toJSON(args, auto_unbox = TRUE)

    handle <- curl::new_handle()
    url <- paste0(url, "/extractRadarGrid")

    curl::handle_setopt(handle, copypostfields = args)
    curl::handle_setheaders(handle, "Content-Type" = "application/json")

    req <- curl::curl_fetch_memory(url, handle = handle)
    if(req$status_code != 200) return(NULL)

    jsonlite::fromJSON(rawToChar(req$content))
}

#' Convert to table extracted radar Cartesian data.
#'
#' Convert to table extracted radar Cartesian data.
#' 
#' @param x A named list, output from \code{extractRadarCart}.
#' 
#' @return A data.frame
#' 
#' @export

gridExtractedTable <- function(x){
    tab <- expand.grid(dates = x$date, 
                       altitude = x$altitude,
                       id = x$coords$id)
    icrd <- match(tab$id, x$coords$id)
    tab$longitude <- x$coords$longitude[icrd]
    tab$latitude <- x$coords$latitude[icrd]
    tab <- tab[, c(3:5, 1:2)]

    for(nm in names(x$data))
        tab[[nm]] <- c(x$data[[nm]])

    return(tab)
}
