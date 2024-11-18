\name{roofDiff}
\alias{roofDiff}
\title{Roof Edge Detection Statistics}
\description{
   Compute the difference between two one-sided gradient estimates.
}
\usage{
   roofDiff(image, bandwidth, blur)
}
\arguments{
    \item{image}{A square matrix, no missing value allowed.}
    \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
    \item{blur}{If blur = TRUE, besides the conventional 2-D kernel
    function, a univariate kernel function is used to address
    the issue of blur.}
}
\value{
   A matrix where each entry is the maximum of the differences:
   \eqn{|\widehat{f}_{x+} - \widehat{f}_{x-}|} and  
   \eqn{|\widehat{f}_{y+} - \widehat{f}_{y-}|} at each pixel location.
}
\details{
   At each pixel, the second-order derivarives (i.e., \eqn{f''_{xx}},
   \eqn{f''_{xy}} and \eqn{f''_{yy}}) are estimated by
   a local quadratic kernel smoothing procedure. Next, the local
   neighborhood is first divided into two halves along the direction
   perpendicular to (\eqn{\widehat{f}''_{xx}}, \eqn{\widehat{f}''_{xy}}).
   Then the one-sided estimates of \eqn{f'_{x+}} and \eqn{f'_{x-}}
   are obtained	respectively by local linear kernel smoothing. The
   estimates of	\eqn{f'_{y+}} and \eqn{f'_{y-}} are obtained by the
   same procedure except that the neighborhood is divided along the
   direction perpendicular to (\eqn{\widehat{f}''_{xy}},
   \eqn{\widehat{f}''_{yy}}).
}
\seealso{
    \code{\link{roofEdgeParSel}}, \code{\link{roofEdge}}
}
\examples{
diff <- roofDiff(image = peppers, bandwidth = 8)
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) "Blind Image Deblurring Using Jump Regression
Analysis", \emph{Statistica Sinica}, \bold{25}, 879 -- 899,
\doi{10.5705/ss.2014.054}.
}