
#' Delete user upload folder
#'
#' Delete user upload folder.
#' 
#' @param dirUPLOAD path to upload folder
#' 
#' @export

deleteUploadDir <- function(dirUPLOAD){
    tz <- "Africa/Kigali"
    Sys.setenv(TZ = tz)
    timeNow <- Sys.time()

    dirs <- list.dirs(dirUPLOAD)
    dirs <- basename(dirs[-1])
    if(length(dirs) == 0) return(NULL)

    temps <- lapply(strsplit(dirs, '_'), function(x) x[length(x)])
    temps <- do.call(c, temps)
    ix <- temps != ""
    if(!any(ix)) return(NULL)
    dirs <- dirs[ix]
    temps <- temps[ix]
    temps <- strptime(temps, "%Y%m%d%H%M%S", tz = tz)
    ## older than 24 hours
    it <- temps < timeNow - 86400
    if(!any(it)) return(NULL)

    del_dir <- file.path(dirUPLOAD, dirs[it])
    unlink(del_dir, recursive = TRUE)
}
