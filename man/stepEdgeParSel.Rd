\name{stepEdgeParSel}
\alias{stepEdgeParSel}
\title{Parameter Selection in Step Edge Detection}
\description{
  Select the bandwidth and threshold parameters for step edge
  detection.
}
\usage{
stepEdgeParSel(image, bandwidth, thresh, nboot, degree = 1,
blur = FALSE)
}
\arguments{
  \item{image}{A square matrix, no missing value allowed.}
  \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
  \item{thresh}{The threshold value to use in the edge detection
    criterion. Must be a positive value.}
  \item{nboot}{Number of bootstrap samples to use in estimating
  \eqn{d_{KQ}}.} 
  \item{degree}{An integer equal to 0 for local constant kernel
    smoothing or 1 for local linear kernel smoothing. The default
    value is 1.}
  \item{blur}{If blur = TRUE, in addition to a conventional 2-D kernel
    function, a 1-D kernel is used in local smoothing to address
    the issue of blur. The default value is FALSE.}
}
\value{
  A list of the selected bandwidth, the selected threshold value
  and a matrix of \eqn{d_{KQ}} values with each entry
  corresponding to each combination of bandwidth and threshold.
}
\details{
    A jump-preserving local linear kernel smoothing is applied to
    estimate the discontinuous regression surface; Bootstrap
    samples are obtained by drawing with replacement from the
    residuals and the \eqn{d_{KQ}} is computed for the detected
    edges of the original sample and those of the bootstrap samples.
}
\seealso{
    \code{\link{roofDiff}}, \code{\link{stepDiff}}, \code{\link{roofEdge}}
}
\examples{
set.seed(24)
parSel <- stepEdgeParSel(image = sar, bandwidth = 5,
thresh = c(17, 21), nboot = 1)
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
