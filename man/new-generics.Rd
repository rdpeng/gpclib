\name{new-generics}
\alias{new-generics}
\alias{append.poly}
\alias{append.poly-methods}
%\alias{append.poly,gpc.poly,gpc.poly-method}
\alias{get.bbox}
\alias{get.bbox-methods}
%\alias{get.bbox,gpc.poly-method}
\alias{area.poly}
\alias{area.poly-methods}
%\alias{area.poly,gpc.poly-method}
\alias{get.pts}
\alias{get.pts-methods}
%\alias{get.pts,gpc.poly-method}
\alias{scale.poly}
\alias{scale.poly-methods}
%\alias{scale.poly,gpc.poly-method}
%\alias{[-methods}
%\alias{[,gpc.poly-method}
\alias{tristrip}
\alias{tristrip-methods}
\alias{triangulate}
\alias{triangulate-methods}

\title{Generics/Methods for polygon objects}

\description{
  Some generic functions and methods for polygon objects
}

\usage{
append.poly(x, y)
area.poly(object, \dots)
get.pts(object)
get.bbox(x)
scale.poly(x, \dots)
tristrip(x)
triangulate(x)
}

\arguments{
  \item{x,object}{A polygon object}
  \item{y}{A polygon object}
  \item{\dots}{Other arguments passed to methods}
}

\section{Methods}{
  \describe{
    \item{append.poly}{\code{signature(x = "gpc.poly", y =
	"gpc.poly")}:  Combine all contours of two \code{"gpc.poly"}
      objects and return the combined polygon as a \code{"gpc.poly"}
      object.}
    \item{area.poly}{\code{signature(object = "gpc.poly")}:  Compute and
      return the sum of the areas of all contours in a \code{"gpc.poly"}
      object.}
    \item{scale.poly}{\code{signature(x = "gpc.poly")}:  Scale (divide)
      the x and y coordinates of a \code{"gpc.poly"} object by the
      amount \code{xscale} and \code{yscale}, respectively.  Return a
      scaled \code{"gpc.poly"} object.}
    \item{get.pts}{\code{signature(object = "gpc.poly")}:  Return the
      list of x and y coordinates of the vertices of a \code{"gpc.poly"}
      object.}
    \item{get.bbox}{\code{signature(x = "gpc.poly")}:  Return the
      bounding box for a \code{"gpc.poly"} object.}
    \item{tristrip}{\code{signature(x = "gpc.poly")}:  Return a tristrip
      list for a  \code{"gpc.poly"} object.}
    \item{triangulate}{\code{signature(x = "gpc.poly")}:  Return a matrix
      of vertices of a triangulation of a \code{"gpc.poly"} object.}
  }
  
}

\details{
The result of \code{tristrip(x)} is a list of two-column matrices.  Each
matrix is a tristrip, i.e. the rows are vertices of triangles, with
each overlapping triple of rows corresponding to a separate triangle.

The result of \code{triangulate(x)} is a single two-column matrix.  The 
rows are vertices of triangles, taken in non-overlapping triples.
  
}


\author{Roger D. Peng; GPC Library by Alan Murta; tristrip additions by
Duncan Murdoch}


\seealso{\code{"gpc.poly"} class documentation.}
\examples{
holepoly <- read.polyfile(system.file("poly-ex/hole-poly.txt", package =
"gpclib"), nohole = FALSE)
area.poly(holepoly)
stopifnot(area.poly(holepoly) == 8)
}
\keyword{methods}
