\name{gpc.poly.nohole-class}
\docType{class}
\alias{gpc.poly.nohole-class}
\alias{coerce,numeric,gpc.poly.nohole-method}

\title{Class "gpc.poly.nohole"}
\description{
  A class for representing polygons with multiple contours but without
  holes.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \sQuote{new("gpc.poly.nohole", ...) or by calling \code{read.polyfile}}.
}
\section{Slots}{
  \describe{
    \item{pts}{Object of class \dQuote{list}.  See the help for
      \dQuote{gpc.poly} for details.}
  }
}
\section{Extends}{
Class \dQuote{gpc.poly}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "numeric", to = "gpc.poly.nohole")}: ... }
  }
}
\author{Roger D. Peng}
\note{
  This class is identical to \dQuote{"gpc.poly"} and is needed because the
  file formats for polygons without holes is slightly different from the
  file format for polygons with holes.  For a \dQuote{gpc.poly.nohole}
  object, the \verb{hole} flag for each contour is always \code{FALSE}.

  Also, \code{write.polyfile} will write the correct file format,
  depending on whether the object is of class \dQuote{gpc.poly} or
  \dQuote{gpc.poly.nohole}.
}

\seealso{
  \code{\link{gpc.poly-class}}
}
\examples{
## None
}
\keyword{classes}
