#' Extract Computed QPE
#'
#' Extract Computed QPE for single scan and hourly
#' 
#' @param start_time,end_time The start and end time same time zone as \code{time_zone}, 
#'                  single scan format "YYYY-mm-dd HH:MM", hourly format "YYYY-mm-dd HH:00"
#' @param dirOUT Full path to folder to save the extracted QPE
#' @param qpe_dir  Full path to folder containing the QPE netCDF files
#' @param qpe_format  File name format of the netCDF files.
#'                    For QPE from single scan, must have 6 "\%s" and 4 "\%s" for hourly data
#' @param extr_type The extraction support, "points" or "shapefiles"
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param padxy If \code{extr_type} is "points", a vector of the padding to use in number of pixels, in order "lon", "lat".
#'              Default c(0, 0), no padding applied.
#' @param shapefiles If \code{extr_type} is "shapefiles", full path to the shapefiles (with the extension ".shp") to use to extract the data
#' @param shp_attr_name If \code{extr_type} is "shapefiles", attribute to be used to extract the data
#' @param shp_attr_values If \code{extr_type} is "shapefiles", vector of the values for \code{shp_attr_name} to extract.
#'                       Default NULL, extract all polygons.
#' @param sp_average If \code{extr_type}  is "shapefiles", if TRUE the extracted data  are spatially average over  the polygons of the shapefiles
#' 
#' @export

extractQPEData <- function(start_time, end_time,
                           dirOUT,
                           qpe_dir = "",
                           qpe_format = "precip_%s%s%s%s%s%s.nc",
                           extr_type = "points",
                           points = NULL,
                           padxy = c(0, 0),
                           shapefiles = NULL,
                           shp_attr_name = NULL,
                           shp_attr_values = NULL,
                           sp_average = TRUE)
{
    if(extr_type == 'points'){
        names(points) <- c('id', 'x', 'y')
        ## pixel size 0.0045 degree ~ 500 m
        padxy <- padxy * 0.0045
        geomObj <- list(type = "points", points = points, padxy = padxy)
    }else{
        dsn <- dirname(shapefiles)
        layer <- tools::file_path_sans_ext(basename(shapefiles))
        shpf <- rgdal::readOGR(dsn = dsn, layer = layer, verbose = FALSE)

        if(!is.null(shp_attr_values)){
            attr_values <- as.character(shp_attr_values)
        }else{
            attr_values <- as.character(shpf@data[, shp_attr_name])
        }
        nom <- gsub("[^[:alnum:]]", "", trimws(attr_values))
        # nom <- stringi::stri_trans_general(str = nom, id = "Latin-ASCII")
        nom <- iconv(nom, from = "UTF-8", to = 'ASCII//TRANSLIT')

        spavg <- if(sp_average) "average" else "gridded"

        geomObj <- list(type = "polys", shp = shpf, attr_values = attr_values,
                        attr_name = shp_attr_name, names = nom, spavg = spavg)
    }

    start <- strptime(start_time, "%Y-%m-%d %H:%M", tz = "UTC")
    end <- strptime(end_time, "%Y-%m-%d %H:%M", tz = "UTC")

    perS <- gregexpr('%s', qpe_format)
    lenS <- length(perS[[1]])
    if(lenS == 6){
        nc_frmt <- sub("%s%s%s%s%s%s", "%s%s%s%s.+\\\\", qpe_format)
        nc_frmt <- paste0(nc_frmt, "$")
        ncInfo <- c("", nc_frmt)
        timestep <- "minute"
        timerange <- list(start = format(start, "%Y-%m-%d-%H-%M"),
                          end = format(end, "%Y-%m-%d-%H-%M"))
    }else{
        ncInfo <- c("", qpe_format)
        timestep <- "hourly"
        timerange <- list(start = format(start, "%Y-%m-%d-%H"),
                          end = format(end, "%Y-%m-%d-%H"))
    }

    ret <- extractQPEnc(timestep, timerange, qpe_dir, ncInfo, dirOUT, geomObj, tz_local = FALSE)

    return(ret)
}

