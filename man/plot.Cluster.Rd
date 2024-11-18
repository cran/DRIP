\name{plot.Surface_Cluster_Parameters}
\alias{plot.Surface_Cluster_Parameters}
\title{Plot Parameter Selection Results in Surface Estimation}
\description{
  Plot information about a clustering-based surface estimation
  parameter selection object.
}
\usage{
\method{plot}{Surface_Cluster_Parameters}(x, ...)
}
\arguments{
  \item{x}{A clustering-based surface estimation parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Plot some information about a clustering-based surface estimation parameter
    selection object. In particular, it plots the cross-validation (no blur) or modified
    cross-validation (there is blur involved) scores against the specified bandwidth
    values.
}
\value{
    A plot of (modified) cross-validation scores is produced. 
}
\seealso{
    \code{\link{surfaceCluster_bandwidth}}, \code{\link{print.Surface_Cluster_Parameters}},
    \code{\link{summary.Surface_Cluster_Parameters}}
}
\examples{
data(brain)
bandwidth_select <- surfaceCluster_bandwidth(image = brain,
    bandwidths = c(3:4), sig.level = .9995, blur = FALSE)
plot(bandwidth_select)
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
