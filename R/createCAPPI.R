
#' Create CAPPI data.
#'
#' Create CAPPI data for single scan radar polar.
#' 
#' @param url The URL of the server. Ex: "http://192.168.1.10:8080"
#' @param source A partial path of the directory containing the folders dates of the mdv files.
#'              Ex: "radarPolar/ops1/sur" or "radarPolar/ops/sur"
#' @param dirOUT Full path to the directory to save the output
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param fields Vector of the fields to create.
#' @param cappi A named list of parameters to be used to create the CAPPI.
#' @param apply_cmd Logical, apply clutter mitigation decision to \code{fields}.
#' @param pia A named list of of the method and parameters to use to perform an attenuation correction for the reflectivity fields \code{dbz_fields}.
#'            Default \code{NULL}, no attenuation correction performed. 
#' @param dbz_fields A vector of reflectivity fields to correct the attenuation. Must be in \code{fields}. Default \code{NULL}
#' @param filter A named list of the method and parameters to use to filter \code{filter_fields}.
#'               Default \code{NULL}, no filtering applied.
#' @param filter_fields A vector of fields in which the filter will be applied. Must be in \code{fields}. Default \code{NULL}
#' @param time_zone the time zone of \code{start_time}, \code{end_time} and the output CAPPI  data.
#'                 Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' 
#' @export

createCAPPI <- function(url, source, dirOUT,
                        start_time, end_time, fields,
                        cappi = list(method = "composite_altitude",
                                     pars = list(fun = "maximum", min_alt = 1.7, max_alt = 15.)),
                        apply_cmd = FALSE,
                        pia = NULL, dbz_fields = NULL,
                        filter = NULL, filter_fields = NULL,
                        time_zone = "Africa/Kigali")
{
    on.exit(curl::handle_reset(handle))
    handle <- curl::new_handle()
    url <- paste0(url, "/createCAPPI")

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
    cappi_pars <- switch(cappi$method,
                         "composite_altitude" = list(fun = "maximum", min_alt = 1.7, max_alt = 15.),
                         "one_altitude" = list(alt = 4.5),
                         "ppi_ranges" = list(alt = 4.5))

    if(!is.null(cappi$pars)){
        if(cappi$method == "composite_altitude"){
            d_name <- names(cappi_pars)
            c_name <- names(cappi$pars)
            inm <- c_name %in% d_name
            if(any(inm)){
                for(n in c_name[inm])
                    cappi_pars[[n]] <- cappi$pars[[n]]
            }
        }else{
            if(!"alt" %in% names(cappi$pars))
                cappi_pars[["alt"]] <- 4.5
        }
    }

    cappi$pars <- cappi_pars

    #######

    start <- strptime(start_time, "%Y-%m-%d %H:%M", tz = time_zone)
    end <- strptime(end_time, "%Y-%m-%d %H:%M", tz = time_zone)
    seqTime <- seq(start, end, "5 min")
    if(time_zone == "UTC"){
        seqTime <- format(seqTime, "%Y-%m-%d-%H-%M")
    }else{
        seqTime <- time_local2utc_char(seqTime, "%Y-%m-%d-%H-%M")
    }

    #######
    
    for(time in seqTime){
        args <- list(source = source,
                     time = time,
                     cappi = cappi,
                     fields = fields,
                     apply_cmd = apply_cmd,
                     pia = pia,
                     dbz_fields = dbz_fields,
                     filter = filter,
                     filter_fields = filter_fields,
                     time_zone = time_zone)
        args <- jsonlite::toJSON(args, auto_unbox = TRUE)

        curl::handle_setopt(handle, copypostfields = args)
        curl::handle_setheaders(handle, "Content-Type" = "application/json")

        req <- curl::curl_fetch_memory(url, handle = handle)
        if(req$status_code != 200){
            msg <- paste("Error occurred, time:", time, time_zone)
            cat(msg, "\n")
            next
        }

        dat <- jsonlite::fromJSON(rawToChar(req$content))
        if(length(dat) == 0){
            msg <- paste("No data, time:", time, time_zone)
            cat(msg, "\n")
            next
        }

        ######
        for(field in fields){
            don <- t(dat$data[[field]])
            dim(don) <- c(length(dat$lon), length(dat$lat), 1)
            dat$data[[field]] <- don
        }

        ######
        ncfile <- paste0("cappi_", dat$time$format, ".nc")
        ncpath <- file.path(dirOUT, ncfile)

        lon <- ncdf4::ncdim_def("lon", "degrees_east", dat$lon, longname = "Longitude")
        lat <- ncdf4::ncdim_def("lat", "degrees_north", dat$lat, longname = "Latitude")
        time <- ncdf4::ncdim_def("time", dat$time$unit, dat$time$value, unlim = TRUE,
                                 calendar = "standard", longname = "Time")

        data_ncout <- vector(mode = "list", length = length(fields))
        for(jj in seq_along(fields)){
            data_ncout[[jj]] <- ncdf4::ncvar_def(fields[jj], "", list(lon, lat, time), -99,
                                                 prec = 'float', longname = fields[jj],
                                                 compression = 6)
        }

        ncout <- ncdf4::nc_create(ncpath, data_ncout)
        for(jj in seq_along(fields))
            ncdf4::ncvar_put(ncout, data_ncout[[jj]], dat$data[[jj]])

        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "title", "Constant Altitude Plan Position Indicator")
        ncdf4::nc_close(ncout)

        cat(paste("Creating CAPPI, time:", dat$time$format, time_zone, "done."), "\n")
    }
}
