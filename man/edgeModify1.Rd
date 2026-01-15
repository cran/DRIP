\name{edgeModify1}
\alias{edgeModify1}
\title{Type-1 Modification Procedure}
\description{
   Modify detected edges to make them thin.
}
\usage{
edgeModify1(bandwidth, image, edge, plot)
}
\arguments{
    \item{image}{A matrix that represents the image.}
    \item{bandwidth}{A positive integer that specifies the number of
    pixels to use in the local smoothing.}
    \item{edge}{A matrix of 0 and 1 represents detected edge
    pixels. }
    \item{plot}{If plot = TRUE, images of detected edges before
    and after the modification are plotted.}
}
\value{
    A matrix of zeros and ones of the same size as the input image.
}
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
    \code{\link{edgeModify2}}
}
\examples{
edge <- stepEdge(sar, bandwidth = 4, thresh = 20, degree = 0)
out <- edgeModify1(4, sar, edge)
}
\author{Yicheng Kang}
\references{
Qiu, P. and Yandell, B. (1997) "Jump Detection in Regression Surfaces,"
\emph{Journal of Computational and Graphical Statistics} \bold{6(3)},
332-354, \doi{10.1080/10618600.1997.10474746}.
}
