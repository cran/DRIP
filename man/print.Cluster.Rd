\name{print.Surface_Cluster_Parameters}
\alias{print.Surface_Cluster_Parameters}
\title{Print Parameter Selection Results in Surface Estimation}
\description{
  Display information about a clustering-based surface estimation
  parameter selection object.
}
\usage{
\method{print}{Surface_Cluster_Parameters}(x, type = "all", ...)
}
\arguments{
  \item{x}{A clustering-based surface estimation parameter selection object.}
  \item{type}{The type of information to display. The "cv_scores"
      option prints the cross-validation or modified cross-validation scores
      for each bandwidth. The "sigma" option prints the estimated noise level.
      The "phi0" option prints the estimated value of the error density at 0.
      The "mean_std_abs" option prints the estimated mean of absolute error.
      The "all" option prints all the information.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Prints some information about a clustering-based surface estimation
    parameter selection object. In particular, this method prints the cross-
    validation or modified cross-validation scores, the selected bandwidth, the
    estimated noise level, the estimated value of the  error density at 0 and
    the estimated mean of absolute error.
}
\value{
    A display of parameter selection results in clustering-based surface
    estimation.
}
\seealso{
    \code{\link{surfaceCluster_bandwidth}},
    \code{\link{summary.Surface_Cluster_Parameters}},
    \code{\link{plot.Surface_Cluster_Parameters}}
}
\examples{
data(brain)
bandwidth_select <- surfaceCluster_bandwidth(image = brain,
    bandwidths = c(3:4), sig.level = .9995, blur = FALSE)
print(bandwidth_select, type = "cv_scores")
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
