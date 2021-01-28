
#' Extract 3D grid data.
#'
#' Extract 3D grid data over a set of points.
#' 
#' @param grid_coords A named list containing the coordinates of the grid data, with names "lon", "lat", "alt"
#' @param grid_data A named list containing the gridded data, the 3d array data are flatten to a vector. 
#'                  The names of list are the names of the data.
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param levels A vector of the index of altitudes to be extracted in integer, or -1 to extract all available altitude
#' @param padxyz A vector of the padding to use in number of pixels, in order "lon", "lat", "alt"
#' @param fun_sp Character, function to be used for the padding. Options: "mean", "median", "max", "min"
#' 
#' @return A named list \cr
#'         List of matrix of the extracted data over the set of points. \cr
#'         The row represents the set of points and the column for the levels.
#' 
#' @export

extract_3DGridData <- function(grid_coords, grid_data,
                               points, levels, padxyz, fun_sp)
{
    grd_lonlat <- list(lon = grid_coords$lon, lat = grid_coords$lat)
    crdGRD <- defSpatialPixels(grd_lonlat)
    ijGRD <- getNeighboursIndex(crdGRD, points, padxyz[1:2])

    nz <- length(grid_coords$alt)
    iz <- lapply(seq_along(levels), function(j){
            z <- levels[j] + (-padxyz[3]:padxyz[3])
            z <- z[z >= 0]
            z <- z[z < nz]
            z + 1
        })

    fun <- get(fun_sp, mode = "function")
    ij2xtr <- lapply(ijGRD$ij2xtr, function(x) x[!is.na(x)])
    nonZero <- sapply(ij2xtr, length) > 0
    ij2xtr <- ij2xtr[nonZero]

    extdat0 <- lapply(seq_along(grid_data), function(j){
        matrix(NA, length(ij2xtr), length(levels))
    })
    names(extdat0) <- names(grid_data)

    if(length(ij2xtr) == 0) return(extdat0)

    ndim <- sapply(grid_coords, length)
    data <- lapply(grid_data, function(field){
        field[is.nan(field)] <- NA
        field <- array(field, ndim)
        array(field, ndim)
    })

    extdat <- lapply(ij2xtr, function(ij){
        don <- lapply(data, function(dfield){
            mat <- lapply(iz, function(i){
                z <- dfield[, , i, drop = FALSE]
                m <- lapply(seq_along(i), function(j) z[, , j][ij])
                m <- do.call(c, m)
                if(length(m) > 1)
                    m <- fun(m, na.rm = TRUE)
                return(m)
            })

            mat <- do.call(c, mat)
            mat[is.nan(mat) | is.infinite(mat)] <- NA
            mat
        })
    })

    nom <- names(data)
    don <- lapply(nom, function(n){
        x <- lapply(extdat, '[[', n)
        do.call(rbind, x)
    })
    names(don) <- nom

    if(!all(nonZero)){
        don <- lapply(seq_along(don), function(j){
            extdat0[[j]][nonZero, ] <- don[[j]]
            extdat0[[j]]
        })
        names(don) <- nom
    }

    return(don)
}
