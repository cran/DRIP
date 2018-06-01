\name{cv.jpex}
\alias{cv.jpex}
\title{Bandwidth Selection and Noise Level Estimation}
\description{
cv.jpex() selects the leave-one-out cross validation (CV) bandwidth for LLK smoothing and estimates the noise level in the input image. Both the bandwidth parameter and the noise level are required inputs for the blind image deblurring procedure jpex().  
}
\usage{
cv.jpex(image, bandwidth)
}
\arguments{
\item{image}{A blurry image to deblurred}

\item{bandwidth}{A vector of positive integers that specifies the size of the neighborhood for local smoothing.}

}
\value{
\item{LLK}{The estimated surface by local linear kernel (LLK) smoothing, using the CV selected bandwidth.}

\item{sigma}{The estimated noise level, defined as the square root of the mean squared error (MSE) between LLK and the input image.}

\item{cv}{A vector of the same length as that of bandwidth. Each element in the vector is the leave-one-out CV error associated with the corresponding bandwidth parameter.}

\item{bandwidth}{The bandwidth parameters input by user.}

\item{band.min}{The bandwidth parameter that results in the smallest CV error.}
}
\references{
Kang, Y. (2018), ``Consistent Blind Image Deblurring Using Jump-Preserving Extrapolation''.
}
\author{Yicheng Kang}
\seealso{\code{\link{jpex}}}
\examples{
require(DRIP)
data(stopsign)
out = cv.jpex(stopsign, c(2,3))
}

