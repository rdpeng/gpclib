\name{polyfile}
\alias{polyfile}
\alias{read.polyfile}
\alias{write.polyfile}

\title{Read/Write polygon data}
\description{
  Read/Write polygon and contour information from/to a text file.  
}
\usage{
read.polyfile(filename, nohole = TRUE)
write.polyfile(poly, filename = "GPCpoly.txt")
}
\arguments{
  \item{filename}{the name of the file (a character string) from/to which
    data should be read/written.}
  \item{nohole}{Is this a polygon without holes?}
  \item{poly}{an object of class \code{"gpc.poly"}}
}
\details{
  The text file representation of a polygon is of the following format:
 
  <number of contours>\cr
  <number of points in first contour>\cr
  x1  y1\cr
  x2  y2\cr
  ...\cr
  <number of points in second contour>\cr
  x1  y1\cr
  x2  y2\cr
  ...\cr

  For example, a data file for a polygon with 2 contours (a four-sided
  object and a triangle) might look like:

  2\cr
  4\cr
  1.0 1.0\cr
  1.0 2.0\cr
  3.4 3.21\cr
  10 11.2\cr
  3\cr
  21.0 11.2\cr
  22.3 99.2\cr
  4.5 5.4\cr

  The vertices of the polygon can be ordered either clockwise or
  counter-clockwise.

  If a polygon has contours which are holes, then the format is slightly
  different.  Basically, a flag is set to indicate that a particular
  contour is a hole.  The format is

  <number of contours>\cr
  <number of points in first contour>\cr
  <hole flag>\cr
  x1  y1\cr
  x2  y2\cr
  ...\cr
  <number of points in second contour>\cr
  <hole flag>\cr
  x1  y1\cr
  x2  y2\cr
  ...\cr

  The hole flag is either 1 to indicate a hole, or 0 for a regular
  contour.  For example, a four-sided polygon with a triangular hole
  would be written as:

  2\cr
  3\cr
  1\cr
  4.0   4.0\cr
  6.0   5.0\cr
  5.0   6.0\cr
  4\cr
  0\cr
  2.0   1.0\cr
  8.0   2.0\cr
  7.0   9.0\cr
  1.0   7.0\cr
}

\value{
  If \code{nohole} is \code{TRUE} (the default) \code{read.polyfile}
  returns an  object of class \code{"gpc.poly.nohole"}.  This object has
  the hole flag set to \code{FALSE} for all contours.  If \code{nohole} is
  \code{FALSE}, then an object of class \code{"gpc.poly"} is
  returned.
  
  \code{write.polyfile} does not return anything useful.
}
\author{Roger D. Peng}

\seealso{
  \code{\link{gpc.poly-class}}, \code{\link{gpc.poly.nohole-class}}
}

\examples{
## None right now.
}
\keyword{IO}
