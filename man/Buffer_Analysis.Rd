% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Buffer_Analysis.R
\name{Buffer_Analysis}
\alias{Buffer_Analysis}
\title{Patch Buffer Analysis}
\usage{
Buffer_Analysis(X, Y = NA, d, scale = 1, minimum_size = 1,
  pathtopython = NULL)
}
\arguments{
\item{X}{wether (i) the path of the file to analyse, (ii) a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).}

\item{Y}{if \code{X} the path of the file to analyse, \code{Y} the path to the folder containing named images of unique colors corresponding to a specific type of organism potentially present on the analysed image. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).}

\item{d}{the distance around each patch that define the buffer.}

\item{scale}{the length of one pixel edge on the image (default to 1).}

\item{minimum_size}{the minimum size of the patches to consider in the analysis (default to 1).}
}
\value{
A dataframe containing the size and relative buffer composition of each patch considered in the analysis.
}
\description{
This function allows you to measure the relative composition of the buffer of a distribution of patch. WARNING: You need the Python library GDAL to use this function. Go to: https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/ for more information on how to proceed.
}
\examples{
#working on it
}
\keyword{Mosaic,}
\keyword{analysis,}
\keyword{buffer,}
\keyword{cover,}
\keyword{distributions.}
\keyword{image}
\keyword{patch}
\keyword{size}
\keyword{size,}
\keyword{species,}
