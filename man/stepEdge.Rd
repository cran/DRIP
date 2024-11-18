\name{stepEdge}
\alias{stepEdge}
\title{Step Edge Detector}
\description{
  Detect step edges in an image.
}
\usage{
stepEdge(image, bandwidth, thresh, degree = 1, blur = FALSE,
plot = FALSE)
}
\arguments{
  \item{image}{A square matrix, no missing value allowed.}
  \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
  \item{thresh}{The threshold value to use in the edge detection
    criterion. Must be a positive value.}
  \item{degree}{An integer equal to 0 for local constant kernel
    smoothing or 1 for local linear kernel smoothing. The default
    value is 1.}
  \item{blur}{If blur = TRUE, in addition to a conventional 2-D kernel
    function, a 1-D kernel is used in local smoothing to address
    the issue of blur. The default value is FALSE.}
  \item{plot}{If plot = TRUE, an image of the detected edges is
  plotted.}
}
\value{
  A matrix of zeros and ones. Ones represent the detected edge
  pixels and zeros represent the non-edge pixels.
}
\details{
    At each pixel, the gradient is estimated by a local	kernel
    smoothing procedure. Next, the local neighborhood is divided
    into two halves along the direction perpendicular to
    (\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the
    one-sided local kernel estimates are obtained in the two half
    neighborhoods respectively. The pixel is flagged as	a step
    edge pixel if \eqn{|\widehat{f}_+ - \widehat{f}_-|>u}, where
    \eqn{u} is the specified threshold value.
}
\seealso{
    \code{\link{roofDiff}}, \code{\link{stepDiff}}, \code{\link{roofEdge}}
}
\examples{
data(sar) 
edges <- stepEdge(image = sar, bandwidth = 4, degree = 0,
thresh = 16)
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
