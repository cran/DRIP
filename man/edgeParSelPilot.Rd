\name{edgeParSelPilot}
\alias{edgeParSelPilot}
\title{Pilot help for parameter selection in edge detection}
\description{
Shows the image of edge detection statistics and upper percentiles for each bandwidth.
}
\usage{
edgeParSelPilot(image, edgeType = "step", degree = 0, blur = FALSE,
                bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95))
}
\arguments{
  \item{image}{The observed image in the form of a matrix.}
  \item{edgeType}{A character: either "step" or "roof" for step edge detection and roof edge detection, respectively.}
  \item{degree}{A numerical value of either 0 or 1. It is the degree of local polynomial smoothing: 0 for local constant kernel smoothing and 1 for local linear kernel smoothing. This argument is not used if edgeType is "roof".}
  \item{blur}{TRUE or FALSE, indicating whether the observed image has blur involved. It affects the kernel function used in the local polynomial smoothing.}
  \item{bandwidth}{A numerical vector specifying the radius of local smoothing in terms of the number of pixels. Only positive integers should be included in the vector.}
  \item{probs}{A numerical vector specifying the probabilities of the edge detection statistics quantiles. These quantiles provide pilot values for the threshold in the edge detection. Usually the upper quantiles are used, as edges are often sparse in an image.}
}
\value{
Returns a matrix with each row and each column corresponding to the bandwidth value and probability. Each cell of the matrix is the quantile of the edge detection statistics with the corresponding bandwidth choice and probability. The function also plots the image of the edge detection statistics.
}
\details{
The bandwidth choice is guided by the visual impression of the images. Note that the edge detection statistics depend the bandwidth, not on the probabilities. After a bandwidth choice is deemed reasonable, the quantiles of the edge detection statistics offers guidance on the threshold values for edge detection.
}
\author{Yicheng Kang}
\references{
Kang, Y. and Qiu, P. (2014) "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}(4), 539 -- 550,
\doi{10.1080/00401706.2013.844732}.
}
\seealso{
\code{\link{roofDiff}}, \code{\link{stepDiff}}
}
\examples{
matout <- edgeParSelPilot(sar, edgeType = "step", degree = 0,
  blur = FALSE, bandwidth = c(5), probs = c(0.8, 0.9))
}
