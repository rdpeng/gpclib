.onLoad <- function(lib, pkg) {
    if(!require(methods, quietly = TRUE))
        stop("'methods' package is required for 'gpclib'")
}

.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- sprintf("General Polygon Clipper Library for R (version %s)\n\tType 'class ? gpc.poly' for help\n", as.character(ver))
    packageStartupMessage(msg)
}
