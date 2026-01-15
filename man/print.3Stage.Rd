\name{print.Three_Stage_Parameters}
\alias{print.Three_Stage_Parameters}
\title{Print Parameter Selection Results in 3-Stage Estimation}
\description{
  Display information about a three-stage estimation bandwidth selection object.
}
\usage{
\method{print}{Three_Stage_Parameters}(x, type = "all", ...)
}
\arguments{
  \item{x}{A three-stage estimation bandwidth selection object.}
  \item{type}{The type of information to display. The "cv_scores"
      option prints the cross-validation or modified cross-validation
      scores for each bandwidth. The "bandwidth" option prints the
      selected the bandwidth. The "all" option prints all the information.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Prints some information about a three-stage estimation bandwidth
    selection object. In particular, this method prints the
    cross-validation or modified cross-validation scores and the
    selected bandwidth.
}
\value{
    Display (modified) cross-validation scores and the selected bandwidth.
}
\seealso{
    \code{\link{restore3StageParSel}}, \code{\link{summary.Three_Stage_Parameters}},
    \code{\link{plot.Three_Stage_Parameters}}
}
\examples{
\dontrun{
step.edges <- stepEdge(peppers, bandwidth = 9, thresh = 17)
roof.edges <- roofEdge(peppers, bandwidth = 6, thresh = 3000,
    edge1 = step.edges)
set.seed(24)
# Time consuming
parSel <- restore3StageParSel(image = peppers, edge1 = step.edges,
    edge2 = roof.edges, bandwidth = 3:4, nboot = 1, blur = TRUE)
print(parSel, type = "all")
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}
