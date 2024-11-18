\name{summary.Three_Stage_Parameters}
\alias{summary.Three_Stage_Parameters}
\title{Summarize Parameter Selection Results in 3-Stage Estimation}
\description{
  Summarize and display key information about a 3-stage estimation bandwidth selection object.
}
\usage{
\method{summary}{Three_Stage_Parameters}(object, ...)
}
\arguments{
  \item{object}{A three-stage estimation bandwidth selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Summarize some information about a three-stage estimation bandwidth
    selection object. In particular, this method prints the
    cross-validation or modified cross-validation scores and the
    selected bandwidth.
}
\value{
    A brief display of bandwidth selection results in three-stage estimation. 
}
\seealso{
    \code{\link{threeStageParSel}}, \code{\link{print.Three_Stage_Parameters}},
    \code{\link{plot.Three_Stage_Parameters}}
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
summary(parSel)
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump 
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}
