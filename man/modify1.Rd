\name{modify1}
\alias{modify1}
\keyword{first-type modification}
\title{Edge detection, post processing}
\description{Modify detected edges to make them thin. }
\usage{modify1(bandwidth, image, edge, plot)}
\arguments{
  \item{image}{A matrix that represents the image.}
  \item{bandwidth}{A positive integer to specify the number of
    pixels used in the local smoothing.}
  \item{edge}{A matrix of 0 and 1 represents detected edge 
    pixels. }
  \item{plot}{If plot = TRUE, images of detected edges before
  the modification and after the modification are plotted.}
}
\value{Returns a matrix of zeros and ones of the same size as
  edge.}
\details{
	A local-smoothing based edge detection algorithm may flag deceptive edge
	pixel candidates. One kind of such candidates consists of those close
	to the real edges. They occur due to the nature of local
	smoothing. That is, if the point \eqn{(x_i, y_j)} is flagged,
	then its neighboring pixels will be flagged with high probability.
	This kind of deceptive candidates can make the detected edges
	thick. This modification procedure makes the detected edges
	thin.
}
\seealso{
	\code{\link{modify2}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a 
          # standard test image in statistics literature.
edge = stepEdgeLCK(sar, 4, 20)
out = modify1(4, sar, edge)
}

\references{
Qiu, P. and Yandell, B., "Jump detection in regression surfaces,"
\emph{Journal of Computational and Graphical Statistics} \bold{6(3)}, 1997,
332-354.
}