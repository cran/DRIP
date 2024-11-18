\name{jpex}
\alias{jpex}
\title{Blind Image Deblurring}
\description{
    Take in any square matrix (noisy blurry image) and deblur it.
}
\usage{
jpex(image, bandwidth, alpha, sigma)
}
\arguments{
    \item{image}{A square matrix representing a blurry image.}
    \item{bandwidth}{A positive integer that specifies the size of
    the neighborhood for local smoothing.}
    \item{alpha}{A numeric between 0 and 1. This is the significance
    level for the Chi-square hypothesis test. The null hypothesis
    is that a given pixel is in a continuity region and not affected
    by the blur.}
    \item{sigma}{A positive numeric value for the noise level in
    the blurred image. It is used in the Chi-square test.}
}
\value{
    \item{deblurred}{A square matrix representing the deblurred
    image.}
    \item{edge}{A square matrix, the element of which is the value
    of the Chi-square test statistic at a pixel location. One can
    classify a given pixel as a blurry pixel if
    \code{edge[i, j] > qchisq(1 - alpha, 2)}.}
}
\references{
Kang, Y. (2020) ``Consistent Blind Image Deblurring Using
Jump-Preserving Extrapolation'', \emph{Journal of Computational and
 Graphical Statistics}, \bold{29}(2), 372 -- 382,
 \doi{10.1080/10618600.2019.1665536}.
}
\author{Yicheng Kang}
\seealso{\code{\link{cv.jpex}}}
\examples{
out <- jpex(image = stopsign, bandwidth = as.integer(2), sigma =
0.00623, alpha = 0.001)
}

