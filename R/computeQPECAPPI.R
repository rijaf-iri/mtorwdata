
#' Compute QPE from CAPPI.
#'
#' Compute QPE from CAPPI for single scan radar polar data.
#' 
#' @param url The URL of the server. Ex: "http://192.168.1.10:8080"
#' @param dirOUT Full path to the directory to save the output
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, format "YYYY-mm-dd HH:MM"
#' @param cappi A named list of parameters to be used to create the CAPPI.
#' @param qpe A named list of parameters to be used to compute the precipitation rate.
#' @param dbz_thres A named list of the minimum and maximum reflectivity threshold.
#'                  Default \code{dbz_thres = list(min = 20, max = 65)}
#' @param apply_cmd Logical, apply clutter mitigation decision to the fields to use to compute the QPE.
#'                  Default \code{TRUE}, applying CMD.
#' @param pia A named list of parameters to be used to correct the attenuation.
#'            Default \code{NULL}, no correction performed. 
#' @param filter A named list of parameters to be used to filter the data.
#'               Default \code{NULL}, no filtering applied.
#' @param time_zone the time zone of \code{start_time}, \code{end_time} and the output QPE.
#'                 Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"
#' 
#' @section CAPPI parameters:
#' cappi
#' 
#' @section QPE parameters:
#' qpe
#'
#' @section Attenuation correction parameters:
#' pia
#' 
#' @section Filter parameters:
#' filter
#' 
#' @return A netCDF files containing the precipitation rate and accumulation
#'         for a single radar scan saved under the folder \code{dirOUT}
#'         with a file name format "precip_YYYYmmddHHMMSS.nc", same time zone as \code{time_zone}.
#' 
#' @export

computeQPECAPPI <- function(url, dirOUT,
                            start_time, end_time, 
                            cappi = list(method = "composite_altitude",
                                         pars = list(fun = "maximum", min_alt = 1.7, max_alt = 15.)),
                            qpe = list(method = "RATE_Z",
                                       pars = list(alpha = 300, beta = 1.4, invCoef = FALSE)),
                            dbz_thres = list(min = 20, max = 65),
                            apply_cmd = TRUE,
                            pia = NULL,
                            filter = NULL,
                            time_zone = "Africa/Kigali")
{
    on.exit(curl::handle_reset(handle))
    handle <- curl::new_handle()
    url <- paste0(url, "/computeQPECAPPI")

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
    qpe_pars <- switch(qpe$method,
                       "RATE_Z" = list(alpha = 300, beta = 1.4, invCoef = FALSE),
                       "RATE_Z_ZDR" = list(alpha = 0.00786, beta_zh = 0.967, beta_zdr = -4.98),
                       "RATE_KDP" = list(alpha = 53.3, beta = 0.669),
                       "RATE_KDP_ZDR" = list(alpha = 192, beta_kdp = 0.946, beta_zdr = -3.45),
                       "RATE_ZPOLY" = NULL)

    if(!is.null(qpe$pars)){
        if(qpe$method == "RATE_ZPOLY"){
            qpe_pars <- NULL
        }else{
            d_name <- names(qpe_pars)
            q_name <- names(qpe$pars)
            inm <- q_name %in% d_name
            if(any(inm)){
                for(n in q_name[inm])
                    qpe_pars[[n]] <- qpe$pars[[n]]
            }
        }
    }

    qpe$pars <- qpe_pars

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
        args <- list(time = time,
                     cappi = cappi,
                     qpe = qpe,
                     dbz_thres = dbz_thres,
                     pia = pia,
                     filter = filter,
                     apply_cmd = apply_cmd,
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
        rate = t(dat$qpe$rate$data)
        rate[rate < 0] <- 0
        dim(rate) <- c(length(dat$lon), length(dat$lat), 1)
        precip = t(dat$qpe$precip$data)
        precip[precip < 0] <- 0
        dim(precip) <- c(length(dat$lon), length(dat$lat), 1)

        ######
        ncfile <- paste0("precip_", dat$time$format, ".nc")
        ncpath <- file.path(dirOUT, ncfile)

        lon <- ncdf4::ncdim_def("lon", "degrees_east", dat$lon, longname = "Longitude")
        lat <- ncdf4::ncdim_def("lat", "degrees_north", dat$lat, longname = "Latitude")
        time <- ncdf4::ncdim_def("time", dat$time$unit, dat$time$value, unlim = TRUE,
                                 calendar = "standard", longname = "Time")
        rate_ncout <- ncdf4::ncvar_def(dat$qpe$rate$name, dat$qpe$rate$unit, list(lon, lat, time), -99,
                                       prec = 'float', longname = dat$qpe$rate$long_name, compression = 6)
        precip_ncout <- ncdf4::ncvar_def(dat$qpe$precip$name, dat$qpe$precip$unit, list(lon, lat, time), -99,
                                         prec = 'float', longname = dat$qpe$precip$long_name, compression = 6)
        data_ncout <- list(rate_ncout, precip_ncout)

        ncout <- ncdf4::nc_create(ncpath, data_ncout)
        ncdf4::ncvar_put(ncout, data_ncout[[1]], rate)
        ncdf4::ncvar_put(ncout, data_ncout[[2]], precip)

        ncdf4::ncatt_put(ncout, "lon", "axis", "X")
        ncdf4::ncatt_put(ncout, "lat", "axis", "Y")
        ncdf4::ncatt_put(ncout, "time", "axis", "T")
        ncdf4::ncatt_put(ncout, 0, "title", "Quantitative Precipitation Estimation")
        ncdf4::nc_close(ncout)

        cat(paste("Computing QPE, time:", dat$time$format, time_zone, "done."), "\n")
    }
}
