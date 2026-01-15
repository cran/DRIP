\name{edgeModify2}
\alias{edgeModify2}
\title{Type-2 Modification Procedure}
\description{
   Delete deceptive edge pixels that are scattered in the design space.
}
\usage{edgeModify2(bandwidth, edge, plot)}
\arguments{
    \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
    \item{edge}{A matrix of 0 and 1 representing the detected edge
    pixels.}
    \item{plot}{If plot = TRUE, images of the detected edges before
    and after the modification are plotted.}
}
\value{
    A matrix of zeros and ones of the same size as the input image.
}
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
    \code{\link{edgeModify1}}
}
\examples{
edge <- stepEdge(sar, bandwidth = 4, thresh = 20, degree = 0)
out <- edgeModify2(4, edge)
}
\author{Yicheng Kang}
\references{
Qiu, P. and Yandell, B. (1997) "Jump Detection in Regression Surfaces,"
\emph{Journal of Computational and Graphical Statistics} \bold{6(3)},
332-354, \doi{10.1080/10618600.1997.10474746}.
}
