\name{diffLC2K}
\alias{diffLC2K}
\keyword{local constant two-kernel estimator}
\title{local constant kernel difference}
\description{Compute difference between two one-sided LC2K
  estimators along the gradient direction.}
\usage{diffLC2K(image, bandwidth, plot)}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{A positive integer to specify the number of
    pixels used in the local smoothing.}
  \item{plot}{If plot = TRUE, an image of the difference at
    each pixel is plotted.}
}
\value{Returns a matrix of the estimated difference, \eqn{|\widehat{f}_+ - \widehat{f}_-|},
	       at each pixel.}
\details{
	At each pixel, the gradient is estimated by a local linear
	kernel smoothing procedure. Next, the local neighborhood is
	divided into two halves along the direction perpendicular to
	(\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the one-
	sided deblurring local constant kernel (LC2K) estimates are obtained in the
	two half neighborhoods respectively.
}
\seealso{
	\code{\link{diffLCK}}, \code{\link{diffLLK}}, \code{\link{diffLL2K}},
	\code{\link{stepEdgeLC2K}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a 
          # standard test image in statistics literature.
diff <- diffLC2K(image = sar, bandwidth = 4)
}

\references{
Kang, Y., and Qiu, P., "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}, 2014, 539-550.
}
