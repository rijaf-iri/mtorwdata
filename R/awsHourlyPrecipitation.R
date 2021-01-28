#' Get hourly precipitation data.
#'
#' Get hourly precipitation data for all AWS.
#' 
#' @param start_time start time, format "YYYY-mm-dd HH:00"
#' @param end_time end time, format "YYYY-mm-dd HH:00"
#' @param url the URL of the server. Ex: "http://192.168.1.10:8080"
#' 
#' @return A named list of coordinates (list name "coords"), dates (list name "date") and the data (list name "data")
#'         \itemize{
#'            \item{\strong{coordinates}: }{a data.frame with column names "id", "longitude" and "latitude"}
#'            \item{\strong{dates}: }{a vector of dates in the format "YYYYmmddHH" in local time}
#'            \item{\strong{data}: }{a matrix with row number equals to the length of dates and column number equals to the length of coordinates}
#'          }
#' 
#' @export

awsGetHourlyPrecip <- function(start_time, end_time, url){
    on.exit(curl::handle_reset(handle))

    data <- list(start_time = start_time, end_time = end_time)
    data <- jsonlite::toJSON(data, auto_unbox = FALSE)

    handle <- curl::new_handle()
    curl::handle_setopt(handle, copypostfields = data)
    curl::handle_setheaders(handle, "Content-Type" = "application/json")

    url <- paste0(url, "/awsGetHourlyPrecip")
    req <- curl::curl_fetch_memory(url, handle = handle)

    out <- jsonlite::fromJSON(rawToChar(req$content))
    return(out)
}

#' Write hourly precipitation data to CDT format.
#'
#' Write hourly precipitation data to CDT format.
#' 
#' @param aws_data Hourly precipitation data output from \code{awsGetHourlyPrecip}
#' @param csvfile Full path to the CSV file to save the CDT format data
#' 
#' @export

awsWriteHourlyPrecip2CDT <- function(aws_data, csvfile){
    dhead <- cbind(c("ID", "LON", "LAT"), t(aws_data$coords))
    don <- cbind(aws_data$date, aws_data$data)
    don <- cbind(t(dhead), t(don))
    don <- t(don)

    write.table(don, csvfile, sep = ",", na = "-99", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
}
