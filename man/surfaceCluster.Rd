\name{surfaceCluster}
\alias{surfaceCluster}
\keyword{local pixel clustering}
\title{Denoising, deblurring and edge-preserving}
\description{ Estimate surface using local pixel clustering 
  and kernel smoothing. Bandwidth is specified by user.}
\usage{surfaceCluster(image, bandwidth, sig.level, sigma, phi0, mean_std_abs, cw=3,
blur = FALSE, plot = FALSE)}
\arguments{
  \item{image}{A square matrix object of size n by n, no missing
    value allowed.}
  \item{bandwidth}{A positive integer that specifies the number of
    pixels used in the local smoothing.}
  \item{sig.level}{Specifies the significance level of the
    hypothesis test deciding to cluster pixels or not.}
  \item{sigma}{Specifies the noise level (i.e., standard deviation of the error distribution). It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated sigma.}
  \item{phi0}{Specifies the density of the standardized error distribution at 0. It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated value.}
  \item{mean_std_abs}{Specifies the mean of absolute value of the standardized error. It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated value.}
  \item{cw}{Specifies the center weight for the center weighted
    median filter. It must be a positive integer.}
  \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
  \item{plot}{ If plot = TRUE, the image of the fitted surface is
    plotted}
}
\value{Returns a list. 'estImg' is the restored image. 'sigma' is the estimated standard deviation of the random error. It is the input value if specified by the user. 'phi0' is the estimated density of the error distribution at 0. It is the input value if specified by the user. 'mean_std_abs' is the estimated absolute mean of the error distribution. It is the input value if specified by the user.}
\seealso{
	\code{\link{JPLLK_surface}}, \code{\link{threeStage}}
}
\examples{
data(brain) 
fit <- surfaceCluster(image=brain, bandwidth=4, sig.level=.9995, cw=3, blur=TRUE)
}
\references{
Kang, Y., Mukherjee, P.S., and Qiu, P. (2017), "Efficient Blind Image
Deblurring Using Nonparametric Regression and Local Pixel 
Clustering", Technometrics, DOI: 10.1080/00401706.2017.1415975.
}
