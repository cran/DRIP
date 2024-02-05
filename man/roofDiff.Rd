\name{roofDiff}
\alias{roofDiff}
\keyword{
	second-order difference estimator
}
\title{
	roof/valley edge detection
}
\description{
	Compute difference between two one-sided gradient
  	estimators.
}
\usage{
	roofDiff(image, bandwidth, blur)
}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{A positive integer to specify the number of
    pixels used in the local smoothing.}
  \item{blur}{If blur = TRUE, besides the conventional 2-D kernel
    function, a univariate kernel function is used to address
    the issue of blur.}
}
\value{
	Returns a matrix where each entry is the maximum of the 
	differences: \eqn{|\widehat{f}_{x+} - \widehat{f}_{x-}|} and  
	\eqn{|\widehat{f}_{y+} - \widehat{f}_{y-}|} at each pixel.
}
\details{
	At each pixel, the second-order derivarives (i.e., \eqn{f''_{xx}},
	\eqn{f''_{xy}}, and \eqn{f''_{yy}}) are estimated by
	a local quadratic kernel smoothing procedure. Next, the local
	neighborhood is first divided into two halves along the direction
	perpendicular to (\eqn{\widehat{f}''_{xx}}, \eqn{\widehat{f}''_{xy}}). Then the
	one-sided estimates of \eqn{f'_{x+}} and \eqn{f'_{x-}} are obtained
	respectively by local linear kernel smoothing. The estimates of
	\eqn{f'_{y+}} and \eqn{f'_{y-}} are obtained by the same procedure
	except that the neighborhood is divided along the direction
	(\eqn{\widehat{f}''_{xy}}, \eqn{\widehat{f}''_{yy}}).
}
\seealso{
	\code{\link{roofEdgeParSel}}, \code{\link{roofEdge}}
}
\examples{
	data(peppers) 
	diff <- roofDiff(image = peppers, bandwidth = 8) # Time consuming
}

\references{
Qiu, P., and Kang, Y. "Blind Image Deblurring Using Jump Regression
Analysis," \emph{Statistica Sinica}, \bold{25}, 2015, 879-899.
}