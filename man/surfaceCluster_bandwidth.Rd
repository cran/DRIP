\name{surfaceCluster_bandwidth}
\alias{surfaceCluster_bandwidth}
\title{Bandwidth Selection for Clustering-Based Surface Estimation}
\description{
    Select the bandwidth parameter for the function
    \code{\link{surfaceCluster}} using cross validation. In the
    cases when there is no blur involved (i.e., denoising only),
    leave-one-out cross validation is used. In the cases when there
    is blur involved, a modified cross validation is used. }
\usage{
surfaceCluster_bandwidth(image, bandwidths, sig.level, sigma,
 phi0, mean_std_abs, relwt = 0.5, cw = 3, blur = FALSE)
}
\arguments{
    \item{image}{A square matrix, no missing value allowed.}
    \item{bandwidths}{An array of positive integers that specifies
    the candidate bandwidth parameters. All the array elements
    must be positive integers because the bandwidth is specified
    in terms of number of pixels.}
    \item{sig.level}{The significance level for the hypothesis
    test deciding whether to cluster pixels or not.}
    \item{sigma}{The noise level (i.e., standard deviation of the
    error distribution). It is used for computing the asymptotic
    threshold for residuals, which are defined to be the
    difference between the local linear kernel smoothing output
    and the center weighted median filter output. If not specified
    by the user, a jump-preserving local linear kernel smoothing
    surface estimation (Qiu 2009) is used to obtain an estimated
    sigma.}
    \item{phi0}{The density of the standardized error distribution
    at 0. It is used for computing the asymptotic threshold for
    residuals, whose definition is given above. If not specified by
    the user, a jump-preserving local linear kernel smoothing
    surface estimation (Qiu 2009) is used to obtain an estimated
    value.}
    \item{mean_std_abs}{The mean of absolute value of the standardized
    error. It is used for computing the asymptotic threshold for
    residuals, whose definition is given above. If not specified by
    the user, a jump-preserving local linear kernel smoothing
    surface estimation (Qiu 2009) is used to obtain an estimated
    value.}
    \item{relwt}{The relative weight assigned to the cross validation
    score in the continuity region. That is, 1 - relwt is assigned
    to the cross validation score around the step edges. It is used
    only when there is blur involved.}
    \item{cw}{The center weight for the center weighted median filter.
    It must be a positive integer.}
    \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
}
\value{
    A list: 'cv_dataframe', a data frame containing the cross
    validation scores corresponding to each candidate bandwidth,
    'bandwidth_hat', the selected bandwidth, 'sigma', the 
    estimated standard deviation of the random error, 'phi0', the
    estimated density of the error distribution at 0, and
    'mean_std_abs', the estimated absolute mean of the error
    distribution.
}
\seealso{
    \code{\link{JPLLK_surface}}, \code{\link{restore3Stage}}
}
\examples{
data(brain) 
bandwidth_select <- surfaceCluster_bandwidth(image = brain,
bandwidths = c(3:4), sig.level = .9995, blur = FALSE)
}
\author{Yicheng Kang}
\references{
Kang, Y., Mukherjee, P.S. and Qiu, P. (2018) "Efficient Blind Image
Deblurring Using Nonparametric Regression and Local Pixel 
Clustering", \emph{Technometrics}, \bold{60}(4), 522 -- 531,
\doi{10.1080/00401706.2017.1415975}.

Qiu, P. (2009) "Jump-Preserving Surface Reconstruction from Noisy 
Data", \emph{Annals of the Institute of Statistical Mathematics},
\bold{61}, 715 -- 751, \doi{10.1007/s10463-007-0166-9}.
}
