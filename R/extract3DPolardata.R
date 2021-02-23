
#' Extract 3D Polar data.
#'
#' Interpolate and extract 3D polar data over a set of points.
#' 
#' @param grid_coords A named list containing the coordinates of the grid data, with names "x", "y", "z"
#' @param grid_data A named list containing the gridded data, the 3d array data are flatten to a vector. 
#'                  The names of list are the names of the data.
#' @param grid_interp A named list containing the coordinates at which the polar data will be interpolated,
#'                    with names "x", "y", "z", "lon" and "lat"
#' @param points A data frame of the points to extract. Data frame with column names "id", "longitude" and "latitude"
#' @param padxy A vector of the padding to use in number of pixels, in order "x", "y"
#' @param fun_sp Character, function to be used for the padding. Options: "mean", "median", "max", "min"
#' @param maxdist The maximum radius of influence in meters for the interpolation. 
#' 
#' @return A named list \cr
#'         List of matrix of the extracted data over the set of points. \cr
#'         The row represents the set of points and the column for the heights.
#' 
#' @export

extract_3DPolarData <- function(grid_coords, grid_data, grid_interp,
                                points, padxy, fun_sp, maxdist = 1000)
{
    fun_sp <- get(fun_sp, mode = "function")
    # d_xyz <- as.data.frame(grid_coords[c('x', 'y', 'z')])
    d_xyz <- as.data.frame(grid_coords)

    lonlat <- list(lon = grid_interp$lon, lat = grid_interp$lat)
    crdGRD <- defSpatialPixels(lonlat)
    ijGRD <- getNeighboursIndex(crdGRD, points, padxy)

    grd_xy <- list(lon = grid_interp$x, lat = grid_interp$y)
    crdXY<- defSpatialPixels(grd_xy)
    crdXYZ <- lapply(ijGRD$ij2xtr, function(ij){
        xy <- sp::coordinates(crdXY[ij])
        ix <- seq(nrow(xy))
        xz <- expand.grid(ix, grid_interp$z)
        xz <- cbind(xy[xz[, 1], , drop = FALSE], xz[, 2])
        xz <- data.frame(xz)
        names(xz) <- c('x', 'y', 'z')
        rownames(xz) <- NULL
        xz
    })

    intrep_data <- lapply(grid_data, function(data){
        data[is.nan(data)] <- NA
        data <- cbind(d_xyz, v = data)

        out <- lapply(crdXYZ, function(xyz){
            index <- split(seq(nrow(xyz)), as.factor(xyz$z))
            rxy <- apply(xyz[, c('x', 'y'), drop = FALSE], 2, range)
            rxy <- maxdist * c(-1, 1) + rxy
            lx <- data$x >= rxy[1, 1] & data$x <= rxy[2, 1]
            ly <- data$y >= rxy[1, 2] & data$y <= rxy[2, 2]

            datap <- data[lx & ly, , drop = FALSE]
            datap <- datap[!is.na(datap$v), , drop = FALSE]

            if(nrow(datap) == 0)
                return(rep(NA, length(index)))

            int <- gstat::idw(v~1, locations = ~x+y+z,
                              data = datap, newdata = xyz,
                              nmax = 8, maxdist = maxdist,
                              debug.level = 0)
            
            pts_elv <- sapply(index, function(ix){
                mat <- int$var1.pred[ix]
                if(length(mat) > 1)
                    mat <- fun_sp(mat, na.rm = TRUE)
                mat
            })

            pts_elv[is.nan(pts_elv) | is.infinite(pts_elv)] <- NA
            pts_elv
        })

        out <- do.call(cbind, out)
        out[is.na(out)] <- NaN
        out
    })

    intrep_data <- lapply(intrep_data, t)

    return(intrep_data)
}
