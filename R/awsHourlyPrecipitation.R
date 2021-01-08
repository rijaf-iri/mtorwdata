#' Get hourly precipitation data.
#'
#' Get hourly precipitation data for all AWS.
#' 
#' @param start_time start time, format "yyyy-mm-dd HH:00"
#' @param end_time end time, format "yyyy-mm-dd HH:00"
#' @param url the URL of the server. Ex: "http://192.168.1.10:8080"
#' 
#' @return a list object
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
