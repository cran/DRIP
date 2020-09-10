\name{stepEdgeLL2K}
\alias{stepEdgeLL2K}
\keyword{local linear two-kernel step edge detection}
\title{Edge detection, denoising and deblurring}
\description{Detect step edges in an image using piecewise
  local linear kernel smoothing.}
\usage{stepEdgeLL2K(image, bandwidth, thresh, plot)}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{A positive integer to specify the number of
    pixels used in the local smoothing.}
  \item{thresh}{Threshold value used in the edge detection
    criterion. }
  \item{plot}{If plot = TRUE, an image of detected edges is
  plotted.}
}
\value{Returns a matrix of zeros and ones of the same size as
  image. Value one represent edge pixels and value zero represent
  non-edge pixels.}
\details{
	At each pixel, the gradient is estimated by a local linear
	kernel smoothing procedure. Next, the local neighborhood is
	divided into two halves along the direction perpendicular to
	(\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the one-
	sided deblurring local linear kernel (LL2K) estimates are obtained in the
	two half neighborhoods respectively. The pixel is flagged as
	a step edge pixel if \eqn{|\widehat{f}_+ - \widehat{f}_-|>u}, where
	\eqn{u} is a threshold value.
}
\seealso{
	\code{\link{stepEdgeLCK}}, \code{\link{stepEdgeLLK}}, \code{\link{stepEdgeLC2K}},
	\code{\link{diffLL2K}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a 
          # standard test image in statistics literature.
edge = stepEdgeLL2K(image = sar, bandwidth = 6, thresh = 20)
}
\references{
Kang, Y., and Qiu, P., "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}, 2014, 539-550.
}