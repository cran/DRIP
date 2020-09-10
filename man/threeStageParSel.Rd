\name{threeStageParSel}
\alias{threeStageParSel}
\keyword{parameter selection in threeStage}
\title{image denoising/deblurring, bandwidth selection, bootstrap}
\description{Select the bandwidth value for the image restoration method
  implemented in the function threeStage}
\usage{threeStageParSel(image, bandwidth, edge1, edge2, nboot, blur=FALSE)}
\arguments{
  \item{image}{A square matrix object of size n by n, no
    missing value allowed.}
  \item{bandwidth}{Bandwidth values to be chosen from. Each of these values need to
  be an positive integer which specifies the number of pixels used in the local smoothing.}
  \item{edge1}{A matrix of 0 and 1 of the same size as
    image represents detected step edge pixels.}
  \item{edge2}{A matrix of 0 and 1 of the same size as
    image represents detected roof/valley edge pixels.}
  \item{nboot}{Required when blur is TRUE. Unused when blur is FALSE. An positive integer to specify the number of bootstraps to perform. See Qiu and Kang (2015) Statistica Sinica for details.}
  \item{blur}{TRUE if the image contains blur, FALSE otherwise. If TRUE, the hybrid selection method proposed in Qiu and Kang (2015) Statistica Sinica is used. If FALSE, the leave-one-out cross validation is used.} 
}
\value{Returns a list of the selected bandwdith, and a matrix of CV values with each entry corresponding to each choice of bandwdith.}
\examples{
data(peppers) # Peppers image is bundled with the package and it is a
          # standard test image in image processing literature.
# Not run
#step.edges = stepEdgeLLK(peppers, 9, 17) # Step edge detection
#roof.edges = roofEdge(peppers, 6, 3000, edge1=step.edges) # Roof edge detection
#parSel = threeStageParSel(image=peppers, edge1=step.edges,
#edge2=roof.edges, bandwidth=4:5) # Time consuming
}
\references{
Qiu, P., and Kang, Y. "Blind Image Deblurring Using Jump Regression
Analysis," \emph{Statistica Sinica}, \bold{25}, 2015, 879-899.
}