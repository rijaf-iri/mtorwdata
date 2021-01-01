
#' Regrid data
#'
#' Regrid data using a bilinear weights.
#' 
#' @param x,y vectors giving the coordinates of the grid to be interpolated
#' @param z vector of the values to be interpolated (same length as x and y)
#' @param x_g,y_g vectors giving the coordinates of the new grid (same length)
#' @param edge logical, if TRUE the grid outside the range of x and y will be extrapolated using the nearest grid
#' 
#' @return a vector of the regridded values
#' 
#' @export

interp_surface_grid <- function(x, y, z, x_g, y_g, edge = TRUE)
{
    loc <- expand.grid(x = x_g, y = y_g)

    nx <- length(x)
    ny <- length(y)

    rule <- if(edge) 2 else 1
    lx <- stats::approx(x, 1:nx, loc$x, rule = rule)$y
    ly <- stats::approx(y, 1:ny, loc$y, rule = rule)$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    ex <- lx - lx1
    ey <- ly - ly1
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1

    z[is.nan(z)] <- NA
    z <- matrix(z, nx, ny)
    z <- z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
         z[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
         z[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
         z[cbind(lx1 + 1, ly1 + 1)] * ex * ey

    return(z)
}

#' Radar Cartesian cross section interpolation
#'
#' Radar Cartesian cross section interpolation.
#' 
#' @param x,y,z vectors giving the coordinates of the grid to be interpolated
#' @param v vector of the values to be interpolated (same length as x, y and z)
#' @param x_g,y_g,z_g vectors giving the coordinates of the new grid (same length)
#' 
#' @return a vector of the interpolated values
#' 
#' @export

interp3d_xsec <- function(x, y, z, v, x_g, y_g, z_g){
    v[is.nan(v)] <- NA
    data <- data.frame(x = x, y = y, z = z, v = v)
    data <- data[!is.na(data$v), , drop = FALSE]
    if(nrow(data) < 4) return(rep(NA, length(x_g)))

    newdata <- data.frame(x = x_g, y = y_g, z = z_g)
    out <- gstat::idw(v~1, locations = ~x+y+z,
                      data = data, newdata = newdata,
                      nmax = 8, maxdist = 1000, debug.level = 0)
    out <- out$var1.pred
    return(out)
}

#' Regrid fields from pseudo CAPPI
#'
#' Regrid fields from pseudo CAPPI
#' 
#' @param x,y vectors giving the coordinates of the grid to be interpolated
#' @param z list of vectors containing the fields to interpolate (same length as x and y)
#' @param x_g,y_g vectors giving the coordinates of the new grid
#' 
#' @return a list of matrices containing the interpolated fields
#' 
#' @export

interp_ppi_ranges_cappi <- function(x, y, z, x_g, y_g){
    nx <- length(x_g)
    ny <- length(y_g)
    grd <- expand.grid(x = x_g, y = y_g)
    coords <- data.frame(x = x, y = y)

    out <- lapply(z, function(v){
        dat <- coords
        v[is.nan(v)] <- NA
        dat$z <- v
        ina <- is.na(dat$z)
        if(all(ina)) return(matrix(NA, nx, ny))
        dat <- dat[!ina, , drop = FALSE]
        kr <- gstat::krige(z~1, locations = ~x+y, data = dat, newdata = grd,
                           nmax = 1, maxdist = 0.0135, debug.level = 0)
        matrix(kr$var1.pred, nx, ny)
    })

    return(out)
}
