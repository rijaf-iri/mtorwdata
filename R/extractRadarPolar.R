
#' Extract radar polar data.
#'
#' Extract radar polar data over a set of points.
#' 
#' @param url The URL of the server. Ex: "http://192.168.1.10:8080"
#' @param source A partial path of the directory containing the folders dates of the mdv files. \cr
#'              Ex: "radarPolar/derived/sur", "radarPolar/ops1/sur" or "radarPolar/ws/sur"
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param fields Vector of the fields to extract.
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param sweeps A vector of the index of elevation angles to be extracted in integer, or -1 to extract all available elevation angles
#' @param pia A named list of of the method and parameters to use to perform an attenuation correction for the reflectivity fields before extraction.
#'            Default \code{NULL}, no attenuation correction performed. 
#' @param dbz_fields A vector of reflectivity fields to correct the attenuation. Must be in \code{fields}. Default \code{NULL}
#' @param filter A named list of the method and parameters to use to filter the fields before extraction.
#'               Default \code{NULL}, no filtering applied.
#' @param filter_fields A vector of fields in which the filter will be applied. Must be in \code{fields}. Default \code{NULL}
#' @param apply_cmd Logical, apply clutter mitigation decision to the fields. Default \code{FALSE}
#' @param time_zone the time zone of \code{start_time}, \code{end_time} and the output QPE.
#'                 Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' 
#' @return A named list
#'\itemize{
#'  \item{\strong{points}: }{the original points used to extract the data}
#'  \item{\strong{date}: }{list of dates of the extracted data}
#'  \item{\strong{elevation_angle}: }{list of the elevation angles of the extracted data}
#'  \item{\strong{data}: }{list of longitude, latitude, altitude and the fields in the form of 3d array.\cr
#'             Dimension: (length(date) x length(elevation_angle) x length(points))}
#' }
#' 
#' @export

extractRadarPolar <- function(url, source,
                              start_time, end_time,
                              fields, points, sweeps = -1,
                              pia = NULL, dbz_fields = NULL,
                              filter = NULL, filter_fields = NULL,
                              apply_cmd = FALSE,
                              time_zone = "Africa/Kigali")
{
    on.exit(curl::handle_reset(handle))

    #######
    if(!is.null(pia)){
        if(pia$method == "kdp"){
            pia_pars = list(gamma = 0.08)
            if(!is.null(pia$pars)){
                if("gamma" %in% names(pia$pars))
                   pia_pars$gamma <- pia$pars$gamma
            }
        }else{
            pia_pars = list(a_max = 0.000167,
                            a_min = 2.33e-05,
                            n_a = 10,
                            b_max = 0.7,
                            b_min = 0.65,
                            n_b = 6,
                            sector_thr = 10,
                            constraints = "none")

            if(!is.null(pia$pars)){
                p_name = names(pia$pars)
                if("constraints" %in% p_name){
                    pia_pars$constraints <- pia$pars$constraints

                    if(pia$pars$constraints == "dbz"){
                        if("constraint_args_dbz" %in% p_name){
                            pia_pars$constraint_args_dbz <- pia$pars$constraint_args_dbz
                        }else{
                            pia_pars$constraint_args_dbz <- 60
                        }
                    }
                    if(pia$pars$constraints == "pia"){
                        if("constraint_args_pia" %in% p_name){
                            pia_pars$constraint_args_pia <- pia$pars$constraint_args_pia
                        }else{
                            pia_pars$constraint_args_pia <- 20
                        }
                    }
                    if(pia$pars$constraints == "both"){
                        if("constraint_args_dbz" %in% p_name){
                            pia_pars$constraint_args_dbz <- pia$pars$constraint_args_dbz
                        }else{
                            pia_pars$constraint_args_dbz <- 60
                        }
                        if("constraint_args_pia" %in% p_name){
                            pia_pars$constraint_args_pia <- pia$pars$constraint_args_pia
                        }else{
                            pia_pars$constraint_args_pia <- 20
                        }
                    }
                }

                d_name <- c("a_max", "a_min", "n_a", "b_max", "b_min", "n_b", "sector_thr")
                inm <- p_name %in% d_name
                if(any(inm)){
                    for(n in p_name[inm])
                        pia_pars[[n]] <- pia$pars[[n]]
                }
            }
        }

        pia$pars <- pia_pars
    }

    #######
    if(!is.null(filter)){
        filter_pars  <- switch(filter$method, 
                               "median_filter_censor" = list(median_filter_len = 5, minsize_seq = 3,
                                                             censor_field = "RHOHV", censor_thres = 0.8),
                               "median_filter" = list(median_filter_len = 5, minsize_seq = 3),
                               "smooth_trim" = list(window_len = 5, window = "hanning")
                             )

        d_name <- names(filter_pars)
        if(!is.null(filter$pars)){
            f_name <- names(filter$pars)
            inm <- f_name %in% d_name
            if(any(inm)){
                for(n in f_name[inm])
                    filter_pars[[n]] <- filter$pars[[n]]
            }

            if(filter$method == "median_filter_censor"){
                if("censor_field" %in% f_name){
                    if(!"censor_thres" %in% f_name){
                        filter_pars[["censor_thres"]] <- switch(filter$pars$censor_field,
                                                                "RHOHV" = 0.8,
                                                                "NCP" = 0.5,
                                                                "SNR" = 3)
                    }
                }
            }
        }

        filter$pars <- filter_pars
    }

    #######

    args <- list(source = source,
                 points = points,
                 start_time = start_time,
                 end_time = end_time,
                 sweeps = sweeps,
                 fields = fields,
                 pia = pia,
                 dbz_fields = dbz_fields,
                 filter = filter,
                 filter_fields = filter_fields,
                 apply_cmd = apply_cmd,
                 time_zone = time_zone
                )
    args <- jsonlite::toJSON(args, auto_unbox = TRUE)

    handle <- curl::new_handle()
    url <- paste0(url, "/extractRadarPolar")

    curl::handle_setopt(handle, copypostfields = args)
    curl::handle_setheaders(handle, "Content-Type" = "application/json")

    req <- curl::curl_fetch_memory(url, handle = handle)
    if(req$status_code != 200) return(NULL)

    jsonlite::fromJSON(rawToChar(req$content))
}


#' Convert to table extracted radar polar data.
#'
#' Convert to table extracted radar polar data.
#' 
#' @param x A named list, output from \code{extractRadarPolar}.
#' 
#' @return A data.frame
#' 
#' @export

polarExtractedTable <- function(x){
    tab <- expand.grid(dates = x$date, 
                       elevation_angle = x$elevation_angle,
                       points_id = x$coords$id)
    icrd <- match(tab$points_id, x$coords$id)
    tab$points_longitude <- x$coords$longitude[icrd]
    tab$points_latitude <- x$coords$latitude[icrd]
    tab <- tab[, c(3:5, 1:2)]

    for(nm in names(x$data))
        tab[[nm]] <- c(x$data[[nm]])

    return(tab)
}
