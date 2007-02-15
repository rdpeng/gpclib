.onLoad <- function(lib, pkg) {
    if(!require(methods, quietly = TRUE))
        stop("'methods' package is required for 'gpclib'")
}

.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- paste("General Polygon Clipper Library for R (version ",
                 as.character(ver), ")", sep = "")
    writeLines(strwrap(msg))
    cat("\tType 'class ? gpc.poly' for help\n")
}
