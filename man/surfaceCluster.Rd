\name{surfaceCluster}
\alias{surfaceCluster}
\title{Jump-Preserving Surface Estimation Using Pixel Clustering}
\description{
    Estimate surface using local pixel clustering and kernel
    smoothing. The bandwidth parameter is specified by the user.
}
\usage{
surfaceCluster(image, bandwidth, sig.level, sigma, phi0,
    mean_std_abs, cw=3, blur = FALSE, plot = FALSE)
}
\arguments{
   \item{image}{A square matrix, no missing   value allowed.}
   \item{bandwidth}{A positive integer that specifies the number
   of pixels to use in the local smoothing.}
   \item{sig.level}{The significance level for the hypothesis
   test deciding whether to cluster pixels or not.}
   \item{sigma}{The noise level (i.e., standard deviation of the
   error distribution). It is used for computing the asymptotic
   threshold for residuals, which are defined to be the difference
   between the local linear kernel smoothing output and the center
   weighted median filter output. If not specified by the user, a
   jump-preserving local linear kernel smoothing surface estimation
   (Qiu 2009) is used to obtain an estimated sigma.}
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
  \item{cw}{The center weight for the center weighted
    median filter. It must be a positive integer.}
  \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
  \item{plot}{ If plot = TRUE, the image of the fitted surface is
    plotted}
}
\value{
    A list of 'estImg', the restored image,  'sigma', the
    estimated standard deviation of the random error, 'phi0', the
    estimated density of the error distribution at 0, and 'mean_std_abs',
    the estimated absolute mean of the error distribution.
}
\seealso{
    \code{\link{JPLLK_surface}}, \code{\link{threeStage}}
}
\examples{ 
fit <- surfaceCluster(image = brain, bandwidth = 4,
    sig.level = .9995, cw = 3, blur = FALSE)
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
