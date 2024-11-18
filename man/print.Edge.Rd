\name{print.Edge_Parameters}
\alias{print.Edge_Parameters}
\title{Print Parameter Selection Results in Edge Detection}
\description{
  Display information about an edge detection parameter selection object.
}
\usage{
\method{print}{Edge_Parameters}(x, type = "all", ...)
}
\arguments{
  \item{x}{An edge detection parameter selection object.}
  \item{type}{The type of information to display. The "matrix"
      option prints the bootstrap version of \code{dKQ} values for each
      combination of bandwidth and threshold. The "parameters"
      option prints the selected bandwidth and threshold. The
      "all" option prints all the information.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Prints some information about an edge detection parameter
    selection object. In particular, this method prints the \eqn{d_{KQ}}
    values based on bootstrap samples, the selected bandwidth and
    threshold.
}
\value{
    A display of parameter selection results in edge detection.
}
\seealso{
    \code{\link{stepEdgeParSel}}, \code{\link{roofEdgeParSel}}
}
\examples{
set.seed(24)
parSel <- stepEdgeParSel(image = sar, bandwidth = 5,
thresh = c(17, 21), nboot = 1)
print(parSel, type = "all")
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
