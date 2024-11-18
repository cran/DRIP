\name{roofEdgeParSel}
\alias{roofEdgeParSel}
\title{Parameter Selection in Roof Edge Detection}
\description{
    Select bandwidth and threshold value for the roof/valley edge detector
    using bootstrap.
}
\usage{
roofEdgeParSel(image, bandwidth, thresh, nboot, edge1, blur = FALSE)
}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{Positive integers to specify the number of
    pixels used in the local smoothing. These are the bandwidth parameters
    to be chosen from.}
  \item{thresh}{Threshold values to be chosen from.}
  \item{nboot}{Number of bootstrap samples.}
  \item{edge1}{Step edges. The function excludes step edges when detect
  roof/valley edges.}
  \item{blur}{TRUE if the image contains blur, FALSE otherwise.} 
}
\value{Returns a list of the selected bandwdith, the selected threshold value,
and a matrix of \eqn{d_{KQ}} values with each entry corresponding to each combination
of bandwdith and threshold.}
\details{
	If \emph{blur=TRUE}, then a conventional local linear kernel smoothing is
	applied to estimate the blurred surface; Bootstrap samples are obtained by
	drawing with replacement from the residuals and the \eqn{d_{KQ}} is computed
	for the detected edges of the original sample and those of the bootstrap
	samples. If \emph{blur=FALSE}, the procedure is the same as when \emph{blur=TRUE}
	except that a jump-preserving kernel smoothing procedure is used to obtain
	residuals.
}
\seealso{
	\code{\link{roofDiff}}, \code{\link{roofEdge}}
}
\examples{
\dontrun{
step.edges <- stepEdge(peppers, bandwidth = 9, thresh = 17)
set.seed(24)
parSel <- roofEdgeParSel(image = peppers, bandwidth = 5, thresh = 5000,
    nboot = 1, edge1 = step.edges, blur = TRUE) # Time Consuming
}
}
\author{Yicheng Kang}
\references{
Qiu, P. and Kang, Y. (2015) ``Blind Image Deblurring Using Jump Regression
Analysis'', \emph{Statistica Sinica}, \bold{25}, 879-899,
\doi{10.5705/ss.2014.054}
}  