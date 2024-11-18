\name{plot.Three_Stage_Parameters}
\alias{plot.Three_Stage_Parameters}
\title{Plot Parameter Selection Results in 3-Stage Estimation}
\description{
  Plot information about a 3-stage estimation bandwidth selection object.
}
\usage{
\method{plot}{Three_Stage_Parameters}(x, ...)
}
\arguments{
  \item{x}{A three-stage estimation bandwidth selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    It plots the cross-validation (no blur involved) or modified cross-validation
    (there is blur involved) scores against the specified bandwidth values.
}
\value{
    A plot of the (modified) cross-validation scores is produced. 
}
\seealso{
    \code{\link{threeStageParSel}}, \code{\link{print.Three_Stage_Parameters}},
    \code{\link{summary.Three_Stage_Parameters}}
}
\examples{
\dontrun{
step.edges <- stepEdge(peppers, bandwidth = 9, thresh = 17) 
roof.edges <- roofEdge(peppers, bandwidth = 6, thresh = 3000,
    edge1 = step.edges) 
set.seed(24)
# Time consuming
parSel <- threeStageParSel(image = peppers, edge1 = step.edges,
    edge2 = roof.edges, bandwidth = 3:4, nboot = 1, blur = TRUE)
plot(parSel)
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump 
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}
