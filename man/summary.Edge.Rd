\name{summary.Edge_Parameters}
\alias{summary.Edge_Parameters}
\title{Summarize Parameter Selection Results in Edge Detection}
\description{
  Summarize and display some key information about an edge detection parameter selection object.
}
\usage{
\method{summary}{Edge_Parameters}(object, ...)
}
\arguments{
  \item{object}{An edge detection parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Summarize results about an edge detection parameter
    selection object. In particular, it shows the \eqn{d_{KQ}}
    values based on bootstrap samples for combinations of bandwidth and
    threshold values.
}
\value{
    A brief display of parameter selection results in edge detection. 
}
\seealso{
    \code{\link{stepEdgeParSel}}, \code{\link{roofEdgeParSel}},
    \code{\link{print.Edge_Parameters}}
}
\examples{
set.seed(24)
parSel <- stepEdgeParSel(image = sar, bandwidth = 5,
thresh = c(17, 21), nboot = 1)
summary(parSel)
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
