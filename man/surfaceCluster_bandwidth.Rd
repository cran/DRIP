\name{surfaceCluster_bandwidth}
\alias{surfaceCluster_bandwidth}
\keyword{bandwidth selection in surfaceCluster}
\title{Denoising, deblurring, bandwidth selection, and edge-preserving}
\description{ Select the bandwidth parameter for the function surfaceCluster based on cross validation. In the cases when there is no blur involved (i.e., denoising only), leave-one-out cross validation is used. In the cases when there is blur involved, a modified cross validation is used. }
\usage{surfaceCluster_bandwidth(image, bandwidths, sig.level, sigma,
 phi0, mean_std_abs, relwt=0.5, cw=3, blur=FALSE)}
\arguments{
  \item{image}{A square matrix object of size n by n, no missing
    value allowed.}
  \item{bandwidths}{An array of positive integers that specifies the candiate bandwidth parameters. All the array elements must be positive integers because the bandwidth is specified in terms of number of pixels.}
  \item{sig.level}{Specifies the significance level of the
    hypothesis test deciding to cluster pixels or not.}
  \item{sigma}{Specifies the noise level (i.e., standard deviation of the error distribution). It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated sigma.}
  \item{phi0}{Specifies the density of the standardized error distribution at 0. It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated value.}
  \item{mean_std_abs}{Specifies the mean of absolute value of the standardized error. It is used for computing the asymptotic threshold for residuals, which are defined to be the difference between the local linear kernel smoothing output and the center weighted median filter output. If not specified by the user, a jump-preserving local linear kernel smoothing surface estimation (Qiu 2009) is used to obtain an estimated value.}
  \item{relwt}{The relative weight assigned to the cross validation score in the continuity region. That is, 1 - relwt is assigned to the cross validation score around the step edges. It is used only when there is blur involved.}
  \item{cw}{Specifies the center weight for the center weighted median filter. It must be a positive integer.}
  \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
}
\value{Returns a list. 'cv_dataframe' contains the cross validation scores corresponding to each candidate bandwidth. 'bandwidth_hat' is the selected bandwidth. 'sigma' is the estimated standard deviation of the random error. It is the input value if specified by the user. 'phi0' is the estimated density of the error distribution at 0. It is the input value if specified by the user. 'mean_std_abs' is the estimated absolute mean of the error distribution. It is the input value if specified by the user. }
\seealso{
	\code{\link{JPLLK_surface}}, \code{\link{threeStage}}
}
\examples{
data(brain) 
bandwidth_select = surfaceCluster_bandwidth(image=brain,
bandwidths=c(3:4), sig.level=.9995, blur=TRUE)
}
\references{
Kang, Y., Mukherjee, P.S., and Qiu, P. (2017), "Efficient Blind Image
Deblurring Using Nonparametric Regression and Local Pixel 
Clustering", Technometrics, DOI: 10.1080/00401706.2017.1415975.

Qiu, P., "Jump-preserving surface reconstruction from noisy data,"
Annals of the Institute of Statistical Mathematics, 61(3), 2009,
715--751.
}
