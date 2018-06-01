\name{jpex}
\alias{jpex}
\title{Blind Image Deblurring}
\description{
jpex() takes in any square matrix (noisy blurry image) and deblurs it.
}
\usage{
jpex(image, bandwidth, alpha, sigma)
}
\arguments{
\item{image}{An input blurry image to deblurred. The input image is represented a square matrix.}

\item{bandwidth}{A positive integer that specifies the size of the neighborhood for local smoothing.}

\item{alpha}{A numberic between 0 and 1. This is the signifance level for the Chi-square hypothesis test that a given pixel is in a continuity region and not affected by the blur.}

\item{sigma}{A positive numeric for the noise level in the blurred image. It is used in the Chi-square test.}
}
\value{
\item{deblurred}{The deblurred image}

\item{edge}{The square matrix, the element of which is the value of the Chi-square test statistic at the pixel location. One can classify a given pixel as a blurry pixcel if edge[i,j]>qchisq(1-alpha, 2).}
}
\references{
Kang, Y. (2018) ``Consistent Blind Image Deblurring Using Jump-Preserving Extrapolation''.
}
\author{Yicheng Kang}
\seealso{\code{\link{cv.jpex}}}
\examples{
require(DRIP)
data(stopsign)
out = jpex(image=stopsign, bandwidth=as.integer(2), sigma=0.00623, alpha=0.001)
}

