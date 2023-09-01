\name{JPLLK_surface}
\alias{JPLLK_surface}
\keyword{Surface estimation}
\title{Denoising and jump-preserving}
\description{ Estimate surface using piecewise local linear 
  kernel smoothing. Bandwidth is chosen by leave-one-out cross validation.}
\usage{JPLLK_surface(image, bandwidth, plot = FALSE)}
\arguments{
  \item{image}{A square matrix object of size n by n, no missing
    value allowed.}
  \item{bandwidth}{A numeric vector with positive integers, which
    specify the number of pixels used in the local smoothing. The 
    final fitted surface chooses the optimal bandwidth from those
    provided by users.}
  \item{plot}{ If plot = TRUE, the image of the fitted surface is
    plotted}
}
\value{A list of fitted values, residuals, chosen bandwidth and
  estimated sigma.}
\details{
	At each pixel, the gradient is estimated by a local linear
	kernel smoothing procedure. Next, the local neighborhood is
	divided into two halves along the direction perpendicular to
	(\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the one-
	sided local linear kernel (LLK) estimates are obtained in the
	two half neighborhoods respectively. Among these two one-sided
	estimates, the one with smaller weighted mean square error is
	chosen to be the final estimate of the regression surface at the
	pixel.
}
\seealso{
	\code{\link{threeStage}}, \code{\link{surfaceCluster}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a 
          # standard test image in statistics literature.
fit <- JPLLK_surface(image=sar, bandwidth=c(3, 4))
}

\references{
Qiu, P., "Jump-preserving surface reconstruction from noisy data", \emph{Annals
of the Institute of Statistical Mathematics}, \bold{61(3)}, 2009, 715-751.
}
