.onLoad <- function(lib, pkg) {
    if(!require(methods, quietly = TRUE))
        stop("'methods' package is required for 'gpclib'")
}

.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- gettextf("General Polygon Clipper Library for R (version %s)",
                    as.character(ver))
    writeLines(strwrap(msg))
    cat(gettext("\tType 'class ? gpc.poly' for help\n"))
}
