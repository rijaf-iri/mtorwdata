
#' Compute validation statistics.
#'
#' A wrapper to compute continuous, categorical and volumetric statistics for validation.
#' 
#' @param stn_data,grd_data A named list of CDT stations data format, output from \code{readCDTStationFormat}.
#' \itemize{
#' \item{\strong{stn_data}: }{contains the reference observations data used to validate}
#' \item{\strong{grd_data}: }{contains the simulated or estimated data to be validated}
#' }
#' @param validation_type Type of the output statistics to compute. Options are
#' \itemize{
#' \item{"stations": }{The statistics will be computed for each stations}
#' \item{"all": }{The statistics will be computed by using all stations together}
#' \item{"average": }{The statistics will be computed by spatially averaging all stations}
#' }
#' @param categorical_thres A number, giving the threshold to use to create the contingency table for the categorical statistics.
#' @param categorical_fun A character, indicating the operator function to be applied for the categorical statistics. Options: ">=", ">", "<=", "<"
#' @param volumetric_type The type of threshold to use to compute the volumetric statistics. Options: "threshold", "percentile"
#' @param volumetric_pars The thresholds or percentile to use to compute the volumetric statistics. If \code{volumetric_type} is
#' \itemize{
#' \item{"threshold": }{\code{volumetric_pars} specifies the threshold to be used.\cr
#' If one value is provided, all the stations will use this threshold\cr
#' If a vector is provided, it must be the same number of the number of stations in \code{stn_data}
#' }
#' \item{"percentile": }{a value of the percentile to use. The percentile will be computed using \code{stn_data} 
#'  to create the threshold to be used to compute the volumetric statistics}
#' }
#' 
#' @return A list object
#' \itemize{
#' \item{\strong{continuous}: }{the continuous statistics. See \code{\link{stats_continuous}} }
#' \item{\strong{categorical}: }{the categorical statistics. See \code{\link{stats_categorical}} }
#' \item{\strong{volumetric}: }{the volumetric statistics. See \code{\link{stats_volumetric_quantile}} }
#' \item{\strong{categorical_pars}: }{Parameters used to compute the categorical statistics.}
#' \item{\strong{volumetric_pars}: }{Parameters used to compute the volumetric statistics.}
#' }
#' 
#' @export

computeValidationStats <- function(stn_data, grd_data, validation_type = "stations", 
                                    categorical_thres = 1, categorical_fun = ">=",
                                    volumetric_type = "percentile", volumetric_pars = 50)
{
    if(!all(stn_data$id == grd_data$id))
        stop("Stations do not match")

    if(!all(stn_data$dates == grd_data$dates))
        stop("Dates do not match")

    if(validation_type == "stations"){
        stn <- stn_data$data
        grd <- grd_data$data
    }

    if(validation_type == "all"){
        stn <- stn_data$data
        stn <- matrix(stn, ncol = 1)
        grd <- grd_data$data
        grd <- matrix(grd, ncol = 1)
    }

    if(validation_type == "average"){
        stn <- stn_data$data
        stn <- rowMeans(stn, na.rm = TRUE)
        stn[is.nan(stn)] <- NA
        stn <- matrix(stn, ncol = 1)
        grd <- grd_data$data
        grd <- rowMeans(grd, na.rm = TRUE)
        grd[is.nan(grd)] <- NA
        grd <- matrix(grd, ncol = 1)
    }

    stat_cont <- stats_continuous(stn, grd)
    pars_categ <- list(thres = categorical_thres, fun = categorical_fun)
    stat_categ <- stats_categorical(stn, grd, pars_categ)

    if(volumetric_type == "threshold"){
        if(length(volumetric_pars) == 1){
            thres <- rep(volumetric_pars, ncol(stn))
        }else{
            thres <- volumetric_pars
        }
        perc_val <- NULL
    }else{
        naPerc <- colSums(!is.na(stn)) < 10
        perc <- volumetric_pars/100
        thres <- matrixStats::colQuantiles(stn, probs = perc, na.rm = TRUE, type = 8)
        thres[naPerc] <- NA
        perc_val <- volumetric_pars
    }

    stat_vol <- stats_volumetric_quantile(stn, grd, thres)
    pars_vol <- list(type = volumetric_type, thres = thres, percentile = perc_val)

    list(continuous = stat_cont,
         categorical = stat_categ,
         categorical_pars = pars_categ,
         volumetric = stat_vol,
         volumetric_pars = pars_vol
     )
}
