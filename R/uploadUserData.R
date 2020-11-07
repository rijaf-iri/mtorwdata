
#' Read csv file
#'
#' Read csv file uploaded by the user.
#' 
#' @param file_csv full path to the csv file
#'
#' @details
#' The data must contain at least 3 column in the order "Name", "Longitude" and "Latitude",
#' and have a header.
#' 
#' @return A json object
#'
#' @export

readUploadCsvCoords <- function(file_csv){
    don <- read.table(file_csv, header = TRUE, sep = ",")
    don <- don[, 1:3, drop = FALSE]
    names(don) <- c("Name", "Longitude", "Latitude")
    return(convJSON(don))
}

#' Read shapefiles
#'
#' Read shapefiles uploaded by the user.
#' 
#' @param dsn full path to the the folder containing files
#' @param layer the files name without extension
#'
#' @details
#' Require the 4 files (.shp, .shx, .dbf and .prj)
#' 
#' @return A geojson object
#'
#' @export

readUploadShptoGeoJSON <- function(dsn, layer){
    shp <- rgdal::readOGR(dsn = dsn, layer = layer, verbose = FALSE)
    shp <- sp::spTransform(shp, sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

    filegeojs <- file.path(dsn, "polygons.geojson")
    rgdal::writeOGR(shp, dsn = filegeojs, layer = layer, driver = "GeoJSON")
    objlist <- jsonlite::fromJSON(filegeojs)
    unlink(filegeojs)

    out <- list(dsn = dsn, layer = layer, data = objlist)
    return(convJSON(out, pretty = FALSE))
}

