\name{stepEdgeParSelLC2K}
\alias{stepEdgeParSelLC2K}
\keyword{edge detection, thresholding}
\title{edge detection, parameter selection}
\description{Select bandwidth and threshold value for LC2K edge detector
  using bootstrap procedure}
\usage{stepEdgeParSelLC2K(image, bandwidth, thresh, nboot)}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{Positive integers to specify the number of
    pixels used in the local smoothing. These are the bandwidth
    parameters to be chosen from.}
  \item{thresh}{Threshold values to be chosen from.}
  \item{nboot}{Number of bootstrap samples.}
}
\value{Returns a list of the selected bandwdith, the selected threshold value,
	and a matrix of \eqn{d_{KQ}} values with each entry corresponding to each
	combination of bandwdith and threshold.}
\details{
	A jump-preserving local linear kernel smoothing is applied to estimate the
	discontinuous regression surface; Bootstrap samples are obtained by
	drawing with replacement from the residuals and the \eqn{d_{KQ}} is computed
	for the detected edges of the original sample and those of the bootstrap
	samples. 
}
\seealso{
	\code{\link{stepEdgeParSelLCK}}, \code{\link{stepEdgeParSelLLK}},
	\code{\link{stepEdgeParSelLL2K}}, \code{\link{stepEdgeLC2K}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a
          # standard test image in statistics literature.
#parSel = stepEdgeParSelLC2K(image=sar, bandwidth=4:5, thresh=
#19:20, nboot=10) # Time consuming
}
\references{
Kang, Y., and Qiu, P., "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}, 2014, 539-550.
}
  