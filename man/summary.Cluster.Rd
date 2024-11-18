\name{summary.Surface_Cluster_Parameters}
\alias{summary.Surface_Cluster_Parameters}
\title{Summarize Parameter Selection Results in Surface Estimation}
\description{
  Summarize and display some key information about a clustering-based surface estimation
  parameter selection object.
}
\usage{
\method{summary}{Surface_Cluster_Parameters}(object, ...)
}
\arguments{
  \item{object}{A clustering-based surface estimation parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Summarize some information about a clustering-based surface estimation parameter
    selection object. In particular, it displays the selected bandwidth, the estimated
    noise level, the estimated value of the  error density at 0 and the estimated
    mean of absolute error.
}
\value{
    A brief display of parameter selection results in clustering-based estimation. 
}
\seealso{
    \code{\link{surfaceCluster_bandwidth}}, \code{\link{print.Surface_Cluster_Parameters}},
    \code{\link{plot.Surface_Cluster_Parameters}}
}
\examples{
data(brain)
bandwidth_select <- surfaceCluster_bandwidth(image = brain,
    bandwidths = c(3:4), sig.level = .9995, blur = FALSE)
summary(bandwidth_select)
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
