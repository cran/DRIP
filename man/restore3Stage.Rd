\name{restore3Stage}
\alias{restore3Stage}
\title{Three-Stage Denoising and Deblurring}
\description{
   Estimate jump location curves using local principal component
   lines. One-sided kernel smoothing is then used for surface
   estimation. Bandwidth is specified by the user.
}
\usage{restore3Stage(image, bandwidth, step_edge = NULL,
roof_edge = NULL, edge1 = NULL, edge2 = NULL, blur = FALSE,
plot = FALSE)
}
\arguments{
    \item{image}{A square matrix, no missing value allowed.}
    \item{bandwidth}{A positive integer that specifies the number
    of pixels to use in the local smoothing.}
    \item{step_edge}{A matrix of 0 and 1 representing the detected
    step edge pixels}
    \item{roof_edge}{A matrix of 0 and 1 representing the detected
    roof/valley edge pixels}
    \item{edge1}{Deprecated. Use \code{step_edge} instead.}
    \item{edge2}{Deprecated. Use \code{roof_edge} instead.}
    \item{blur}{If blur = TRUE, besides a conventional 2-D kernel
    function, a univariate increasing kernel function is used in
    the local kernel smoothing to address the issue with blur.}
    \item{plot}{ If plot = TRUE, the image of the fitted surface is
    plotted}
}
\value{
   The restored image, which is represented by a matrix.
}
\details{
    At each pixel, if there are step edges detected in the local
    neighborhood, a principal component line is fitted through the
    detected edge pixels to approximate	the step edge locally and
    then the regression surface is estimated by a local constant
    kernel smoothing procedure using only the pixels on one side
    of the principal component line. If there are no step edges
    but roof/valley edges detected in the local neighborhood, the
    same procedure is followed except that the principal component
    line to fitted through the detected roof/valley edge pixels.
    In cases when there is either no step edges or roof/valley
    edges detected in the neighborhood, the regression surface at
    the pixel is estimated by the conventional local linear kernel
    smoothing procedure.
}
\seealso{
    \code{\link{JPLLK_surface}}, \code{\link{surfaceCluster}}
}
\examples{
step.edges <- stepEdge(sar, bandwidth = 4, thresh = 20, degree = 0)
stepEdge1 <- edgeModify2(bandwidth = 4, step.edges)
fit <- restore3Stage(image = sar, bandwidth = 4, step_edge = stepEdge1,
    roof_edge = array(0, rep(ncol(sar), 2)))
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump
Regression Analysis'', \emph{Statistica Sinica}, \bold{25},
879 -- 899, \doi{10.5705/ss.2014.054}.
}
