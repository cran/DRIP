\name{modify2}
\alias{modify2}
\keyword{post processing, edge detection}
\title{Edge detection, post processing}
\description{Delete deceptive edge pixels that are scattered
  in the design space.}
\usage{modify2(bandwidth, edge, plot)}
\arguments{
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
	pixel candidates. One kind of such candidates existis due to the nature of hypothesis
	testing, on which the threshold value of the edge detection criterion is
	based. That is, a point \eqn{(x_i, y_j)} could be flagged as a edge
	pixel with certain probability, even if it is actually not a edge pixel.
	Deceptive candidates of this kind are scattered in the whole design
	space. This modification procedure deletes scattered edge pixel candidates.
}

\seealso{
	\code{\link{modify1}}
}
\examples{
data(sar) # SAR image is bundled with the package and it is a 
          # standard test image in statistics literature.
edge = stepEdgeLCK(sar, 4, 20)
out = modify2(4, edge)
}

\references{
Qiu, P. and Yandell, B., "Jump detection in regression surfaces,"
\emph{Journal of Computational and Graphical Statistics} \bold{6(3)}, 1997,
332-354.
}