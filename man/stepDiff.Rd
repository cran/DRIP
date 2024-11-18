\name{stepDiff}
\alias{stepDiff}
\title{Step Edge Detection Statistics}
\description{
  Compute difference between two one-sided local kernel
  estimates along the gradient direction.
}
\usage{
stepDiff(image, bandwidth, degree = 1, blur = FALSE, plot = FALSE)
}
\arguments{
  \item{image}{A square matrix, no missing value allowed.}
  \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
  \item{degree}{An integer equal to 0 for local constant kernel
    smoothing or 1 for local linear kernel smoothing. The default
    value is 1.}
  \item{blur}{If blur = TRUE, in addition to a conventional 2-D kernel
    function, a 1-D kernel is used in local smoothing to address
    the issue of blur. The default value is FALSE.}
  \item{plot}{If plot = TRUE, an image of the detection statistics at
    each pixel is plotted.}
}
\value{
  A matrix of the estimated difference, \eqn{|\widehat{f}_+ - \widehat{f}_-|},
  at each pixel.}
\details{
	At each pixel, the gradient is estimated by a local linear
	kernel smoothing procedure. Next, the local neighborhood is
	divided into two halves along the direction perpendicular to
	(\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the one-
	sided local kernel estimates are obtained in the two half
	neighborhoods respectively.
}
\seealso{
    \code{\link{roofDiff}}
}
\examples{
diff <- stepDiff(image = sar, bandwidth = 4, degree = 0)
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
