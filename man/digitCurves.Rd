% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digit_curves.R
\name{digitCurves}
\alias{digitCurves}
\title{Mesh Digitization}
\usage{
digitCurves(M, ...)
}
\arguments{
\item{M}{Either a mesh3d object (in this case user should refer to \code{\link{digitCurves.mesh3d}}).}

\item{...}{Additional arguments (all are not optional!) needed for curve digitization.}
}
\value{
A numerical matrix (\code{\link{digitCurves.mesh3d}}) with attributes.
}
\description{
Generic function for mesh digitization. It invokes 2 particular methods depending on the class of the
             1st argument \code{M}. So far method is implemented only for the mesh3d object.
}
\details{
For details, user should refer to the method \code{\link{digitCurve.mesh3d}} (to digitize a single
         \code{mesh3d} object). \cr
             \cr
             \strong{WARNING}: For Mac users, \code{digitCurves} is currently not compatible with the RStudio interface.
                               Please, use the basic R interface instead. A call from RStudio will cause an error
                               and an exit from the function... \cr
}
\seealso{
\code{\link{digitCurves.mesh3d}}.
}
