\name{cv.jpex}
\alias{cv.jpex}
\title{Bandwidth Selection and Noise Level Estimation}
\description{
    Select the leave-one-out cross validation bandwidth for local
    linear kernel smoothing and estimates the noise level in the
    input image. Both the bandwidth parameter and the noise level
    are required inputs for the blind image deblurring procedure
    \code{jpex}.  
}
\usage{
cv.jpex(image, bandwidths, ncpus = 1)
}
\arguments{
    \item{image}{A blurry input image.}
    \item{bandwidths}{A vector of positive integers that specify the
    size  of the neighborhood for local smoothing.}
    \item{ncpus}{The number of CPUs allocated for parallel computing.}
}
\value{
    \item{LLK}{The estimated surface by local linear kernel (LLK)
    smoothing, using the CV selected bandwidth.}
    \item{sigma}{The estimated noise level, defined as the square root
    of the mean squared error (MSE) between LLK and the input image.}
    \item{cv_scores}{A vector of the same length as that of the input bandwidths.
    Each element in the vector is the leave-one-out CV score associated
    with the corresponding bandwidth value.}
    \item{input}{The bandwidth parameters input by user.}
    \item{bandwidth}{The bandwidth parameter that results in the smallest
    CV score.}
}
\seealso{\code{\link{jpex}}}
\examples{
out <- cv.jpex(stopsign, c(2,3))
}
\author{Yicheng Kang}
\references{
Kang, Y. (2020) ``Consistent Blind Image Deblurring Using
Jump-Preserving Extrapolation'', \emph{Journal of Computational and
 Graphical Statistics}, \bold{29}(2), 372 -- 382,
 \doi{10.1080/10618600.2019.1665536}.
}



