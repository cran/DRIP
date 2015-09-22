\name{surfaceCluster}
\alias{surfaceCluster}
\keyword{Surface estimation, local pixel clustering}
\title{Denoising, deblurring and edge-preserving}
\description{ Estimate surface using local pixel clustering 
  and kernel smoothing. Bandwidth is specified by user.}
\usage{surfaceCluster(image, bandwidth, sig.level,
blur = FALSE, plot = FALSE)}
\arguments{
  \item{image}{A square matrix object of size n by n, no missing
    value allowed.}
  \item{bandwidth}{A positive integer that specifies the number of
    pixels used in the local smoothing.}
  \item{sig.level}{Specifies the significance level of the
    hypothesis test deciding to cluster pixels or not.}
  \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
  \item{plot}{ If plot = TRUE, the image of the fitted surface is
    plotted}
}
\value{Returns a list. 'estImg' is the restored image. 'sigma2' is
  the estimated second moment of the random error. 'mmt4' is the
  estimated fourth moment of the random error.}
\seealso{
	\code{\link{JPLLK_surface}}, \code{\link{threeStage}}
}
\examples{
data(brain) 
fit = surfaceCluster(image=brain, bandwidth=4, sig.level=.995,
blur=TRUE)
}
\references{
Kang, Y., Mukherjee, P.S., and Qiu, P., "Efficient Blind Image
Deblurring Using Nonparametric Regression and Local Pixel 
Clustering", \emph{Computational Statistics and Data Analysis},
submitted.
}
