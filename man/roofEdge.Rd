\name{roofEdge}
\alias{roofEdge}
\title{Roof Edge Detector}
\description{
    Detect roof/valley edges in an image using piecewise local
    quadratic kernel smoothing.
}
\usage{
roofEdge(image, bandwidth, thresh, edge1, blur, plot)
}
\arguments{
    \item{image}{A square matrix, no missing value allowed.}
    \item{bandwidth}{A positive integer that specifies the number
    of pixels to use in the local smoothing.}
    \item{thresh}{Threshold value to use in the edge detection
    criterion.}
    \item{edge1}{A square matrix representing the image's step
    edges. The function excludes step edges when detects roof
    edges.}
    \item{blur}{If blur = TRUE, besides the conventional 2-D
    kernel function, a univariate kernel function is used
    in the local smoothing to address the issue of blur.}
    \item{plot}{If plot = TRUE, an image of detected edges is
    plotted.}
}
\value{
    A matrix of zeros and ones of the same size as the input image.
}
\details{
    At each pixel, the second-order derivarives (i.e., \eqn{f''_{xx}},
    \eqn{f''_{xy}}, and \eqn{f''_{yy}}) are estimated by a local
    quadratic kernel smoothing procedure. Next, the local neighborhood
    is first divided into two halves along the direction perpendicular
    to (\eqn{\widehat{f}''_{xx}}, \eqn{\widehat{f}''_{xy}}). Then the
    one-sided estimates of \eqn{f'_{x+}} and \eqn{f'_{x-}} are obtained
    respectively by local linear kernel smoothing. The estimates of
    \eqn{f'_{y+}} and \eqn{f'_{y-}} are obtained by the same procedure
    except that the neighborhood is divided along the direction perpendicular
    to (\eqn{\widehat{f}''_{xy}}, \eqn{\widehat{f}''_{yy}}). The pixel is
    flagged as a roof/valley edge pixel if \eqn{max(|\widehat{f}_{x+} -
    \widehat{f}_{x-}|,  |\widehat{f}_{y+} - \widehat{f}_{y-}|)>} the
    specified threshold and there is no step edge pixels in the neighborhood.
}
\seealso{
    \code{\link{roofEdgeParSel}}, \code{\link{roofDiff}}
}
\examples{
\dontrun{
step.edges <- stepEdge(peppers, bandwidth = 6, thresh = 25, degree = 1)
roof.edges <- roofEdge(image = peppers, bandwidth = 9, thresh = 3000,
    edge1 = step.edges, blur = FALSE, plot = FALSE) # Time consuming
edges <- step.edges + roof.edges     
par(mfrow = c(2, 2))
image(1 - step.edges, col = gray(0:1))
image(1 - roof.edges, col = gray(0:1))
image(1 - edges, col = gray(0:1))
image(peppers, col = gray(c(0:255)/255))
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump 
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}