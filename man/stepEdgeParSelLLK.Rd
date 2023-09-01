\name{stepEdgeParSelLLK}
\alias{stepEdgeParSelLLK}
\keyword{LLK parameter selection}
\title{edge detection, parameter selection}
\description{Select threshold value for LLK edge detector
  using bootstrap procedure}
\usage{stepEdgeParSelLLK(image, bandwidth, thresh, nboot)}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{Positive integers to specify the number of
    pixels used in the local smoothing. These are the bandwidth
    parameters to be chosen from.}
  \item{thresh}{Threshold values to be chosen from.}
  \item{nboot}{Number of bootstrap samples.}
}
\value{Returns a list of the selected bandwdith, the selected
threshold value, and a matrix of \eqn{d_{KQ}} values with each entry
corresponding to each combination of bandwdith and threshold.}
\details{
	A jump-preserving local linear kernel smoothing is applied to estimate the
	discontinuous regression surface; Bootstrap samples are obtained by
	drawing with replacement from the residuals and the \eqn{d_{KQ}} is computed
	for the detected edges of the original sample and those of the bootstrap
	samples. 
}
\seealso{
	\code{\link{stepEdgeParSelLCK}}, \code{\link{stepEdgeParSelLC2K}},
	\code{\link{stepEdgeParSelLL2K}}, \code{\link{stepEdgeLLK}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a
          # standard test image in statistics literature.
set.seed(24)
parSel <- stepEdgeParSelLLK(image = sar, bandwidth = 5, thresh = c(17,21), nboot = 1)
}
\references{
Kang, Y., and Qiu, P., "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}, 2014, 539-550.
}