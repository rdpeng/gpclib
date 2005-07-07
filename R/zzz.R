.onLoad <- function(lib, pkg) {
    if(!require(methods))
        stop(sQuote("methods"), " package is required for ", sQuote("gpclib"))
}

.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- paste("General Polygon Clipper Library for R (version ",
                 as.character(ver), ")", sep = "")
    writeLines(strwrap(msg))
    cat(gettext("\tType", sQuote("class ? gpc.poly"), "for help\n"))
}
