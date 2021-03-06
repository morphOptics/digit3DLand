% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{save.tps}
\alias{save.tps}
\title{Saves 3D Coordinates into a TPS File}
\usage{
save.tps(A, ID, file.name, LMheader = "LM3", IDheader = "ID",
         sdir = getwd(), app = FALSE, over.write = FALSE)
}
\arguments{
\item{A}{A 2D numerical matrix of landmark coordinates.}

\item{ID}{A character value containing the individual identifier.}

\item{file.name}{A character value containing the TPS file name (with extension).}

\item{LMheader}{A character value specifying the keyword to use in the TPS file to indicate the number of
landmarks. \cr
Default: \code{"LM3"}, for 3D landmarks.}

\item{IDheader}{A character value specifying the keyword to use in the TPS file to indicate the individual
identifier. \cr
Default: \code{"ID"}.}

\item{sdir}{A character value indicating the path to the saving directory. \cr
Default: \code{getwd()}.}

\item{app}{A logical value indicating if the written data should be appended in the TPS file (if this file already
exists). \cr
Default: \code{FALSE}.}

\item{over.write}{A logical value indicating if the file to save should overwrite the previous TPS file with the
same name (if it exists). \cr
Default: \code{FALSE}.}
}
\value{
Nothing.
}
\description{
Saves landmark coordinates as generated by \code{\link{digitMesh.mesh3d}} into a TPS file. \cr
             \strong{Note}: this is not a generic function to create TPS file, it's only designed for 3D landmark
                            data.
}
\details{
Given a matrix of landmark coordinates \code{A}, this function saves it into a TPS file whose filename is
         contained into \code{file.name}.
}
\examples{

# Generates first a 10*3 random matrix, just for illustration purpose,
# but such as it could be obtained after a mesh digitization:
M <- matrix(rnorm(30), 10, 3)
# Then, saves it as a TPS file:
save.tps(M, ID = "randomMatrix", file.name = "TPSfile4randomMatrix.tps")


}
\seealso{
\code{\link{read.tps}}.
}
