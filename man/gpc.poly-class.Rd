\name{gpc.poly-class}
\docType{class}
\alias{gpc.poly-class}
\alias{[,gpc.poly-method}
\alias{append.poly,gpc.poly,gpc.poly-method}
\alias{area.poly,gpc.poly-method}
\alias{coerce,matrix,gpc.poly-method}
\alias{coerce,data.frame,gpc.poly-method}
\alias{coerce,gpc.poly,matrix-method}
\alias{coerce,gpc.poly,numeric-method}
\alias{coerce,numeric,gpc.poly-method}
\alias{get.bbox,gpc.poly-method}
\alias{get.pts,gpc.poly-method}
\alias{intersect,gpc.poly,gpc.poly-method}
\alias{plot,gpc.poly-method}
\alias{scale.poly,gpc.poly-method}
\alias{setdiff,gpc.poly,gpc.poly-method}
\alias{show,gpc.poly-method}
\alias{union,gpc.poly,gpc.poly-method}
\alias{tristrip,gpc.poly-method}
\alias{triangulate,gpc.poly-method}

\title{Class "gpc.poly"}
\description{
  A class for representing polygons composed of multiple contours, some
  of which may be holes.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("gpc.poly",
    ...)} or by reading in from a file using \code{read.polyfile}.
}
\section{Slots}{
  \describe{
    \item{\code{pts}:}{Object of class \code{"list"}.  Actually,
      \code{pts} is a list of lists with length equal to the number of
      contours in the \code{"gpc.poly"} object.  Each element of
      \code{pts} is a list of length 3 with names \code{x}, \code{y},
      and \code{hole}.  \code{x} and \code{y} are vectors containing the
      x and y coordinates, respectively, while \code{hole} is a logical
      indicating whether or not the contour is a hole.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "gpc.poly")}: ... }
    \item{append.poly}{\code{signature(x = "gpc.poly", y = "gpc.poly")}: ... }
    \item{area.poly}{\code{signature(object = "gpc.poly")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "gpc.poly")}: ... }
    \item{coerce}{\code{signature(from = "data.frame", to = "gpc.poly")}: ... }
    \item{coerce}{\code{signature(from = "gpc.poly", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "gpc.poly", to = "numeric")}: ... }
    \item{coerce}{\code{signature(from = "numeric", to = "gpc.poly")}: ... }
    \item{get.bbox}{\code{signature(x = "gpc.poly")}: ... }
    \item{get.pts}{\code{signature(object = "gpc.poly")}: ... }
    \item{intersect}{\code{signature(x = "gpc.poly", y = "gpc.poly")}: ... }
    \item{plot}{\code{signature(x = "gpc.poly")}: The argument
      \code{poly.args} can be used to pass a list of additional
      arguments to be passed to the underlying \code{polygon} call.}
    \item{scale.poly}{\code{signature(x = "gpc.poly")}: ... }
    \item{setdiff}{\code{signature(x = "gpc.poly", y = "gpc.poly")}: ... }
    \item{show}{\code{signature(object = "gpc.poly")}: Scale x and y
      coordinates by amount \code{xscale} and \code{yscale}.  By default
      \code{xscale} equals \code{yscale}.}
    \item{union}{\code{signature(x = "gpc.poly", y = "gpc.poly")}: ... }
    \item{tristrip}{\code{signature(x = "gpc.poly")}: ... }
    \item{triangulate}{\code{signature(x = "gpc.poly")}: ... }


  }
}
\author{Roger D. Peng}

\note{
  The class \code{"gpc.poly.nohole"} is identical to
  \code{"gpc.poly"} except the \code{hole} flag for each contour of a
  \code{"gpc.poly.nohole"} object is always \code{FALSE}.
}

\examples{
## Make some random polygons
set.seed(100)
a <- cbind(rnorm(100), rnorm(100))
a <- a[chull(a), ]

## Convert `a' from matrix to "gpc.poly"
a <- as(a, "gpc.poly")

b <- cbind(rnorm(100), rnorm(100))
b <- as(b[chull(b), ], "gpc.poly")

## More complex polygons with an intersection
p1 <- read.polyfile(system.file("poly-ex/ex-poly1.txt", package = "gpclib"))
p2 <- read.polyfile(system.file("poly-ex/ex-poly2.txt", package = "gpclib"))

## Plot both polygons and highlight their intersection in red
plot(append.poly(p1, p2))
plot(intersect(p1, p2), poly.args = list(col = 2), add = TRUE)

## Highlight the difference p1 \ p2 in green
plot(setdiff(p1, p2), poly.args = list(col = 3), add = TRUE)

## Highlight the difference p2 \ p1 in blue
plot(setdiff(p2, p1), poly.args = list(col = 4), add = TRUE)

## Plot the union of the two polygons
plot(union(p1, p2))

## Take the non-intersect portions and create a new polygon
## combining the two contours
p.comb <- append.poly(setdiff(p1, p2), setdiff(p2, p1))
plot(p.comb, poly.args = list(col = 2, border = 0))

## Coerce from a matrix
x <- 
structure(c(0.0934073560027759, 0.192713393476752, 0.410062456627342, 
0.470020818875781, 0.41380985426787, 0.271408743927828, 0.100902151283831, 
0.0465648854961832, 0.63981588032221, 0.772382048331416,
0.753739930955121, 0.637744533947066, 0.455466052934407,
0.335327963176065, 0.399539700805524, 
0.600460299194476), .Dim = c(8, 2))
y <- 
structure(c(0.404441360166551, 0.338861901457321, 0.301387925052047, 
0.404441360166551, 0.531852879944483, 0.60117973629424, 0.625537820957668, 
0.179976985040276, 0.341542002301496, 0.445109321058688,
0.610817031070196, 0.596317606444189, 0.459608745684695,
0.215189873417722), .Dim = c(7, 2))

x1 <- as(x, "gpc.poly")
y1 <- as(y, "gpc.poly")

plot(append.poly(x1, y1))
plot(intersect(x1, y1), poly.args = list(col = 2), add = TRUE)

## Show the triangulation
plot(append.poly(x1, y1))
triangles <- triangulate(append.poly(x1,y1))
for (i in 0:(nrow(triangles)/3 - 1)) 
    polygon(triangles[3*i + 1:3,], col="lightblue")


}
\keyword{classes}
