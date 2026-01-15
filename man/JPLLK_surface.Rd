\name{JPLLK_surface}
\alias{JPLLK_surface}
\title{Jump-Preserving Local Linear Kernel Smoothing}
\description{
    Estimate surface using piecewise local linear kernel smoothing.
    The bandwidth is chosen by leave-one-out cross validation.
}
\usage{
JPLLK_surface(image, bandwidth, plot = FALSE)
}
\arguments{
    \item{image}{A square matrix, no missing value allowed.}
    \item{bandwidth}{A numeric vector of positive integers, which
    specifies the number of pixels used in the local smoothing. The 
    final fitted surface uses the optimal bandwidth chosen from
    those provided by users.}
    \item{plot}{If plot = TRUE, the image of the fitted surface is
    plotted.}
}
\value{
    A list of fitted values, residuals, chosen bandwidth and
    estimated sigma.
}
\details{
    At each pixel, the gradient is estimated by a local linear
    kernel smoothing procedure. Next, the local neighborhood is
    divided into two halves along the direction perpendicular to
    (\eqn{\widehat{f}'_{x}}, \eqn{\widehat{f}'_{y}}). Then the one-
    sided local linear kernel (LLK) estimates are obtained in the
    two half neighborhoods respectively. Among these two one-sided
    estimates, the one with smaller weighted mean square error is
    chosen to be the final estimate of the regression surface at the
    pixel.
}
\seealso{
    \code{\link{restore3Stage}}, \code{\link{surfaceCluster}}
}
\examples{
fit <- JPLLK_surface(image = sar, bandwidth = c(3, 4))
}
\author{Yicheng Kang}
\references{
Qiu, P. (2009) "Jump-Preserving Surface Reconstruction from Noisy Data",
\emph{Annals of the Institute of Statistical Mathematics}, \bold{61}(3),
715 -- 751, \doi{10.1007/s10463-007-0166-9}.
}
