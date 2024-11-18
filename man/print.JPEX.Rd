\name{print.JPEX_Parameters}
\alias{print.JPEX_Parameters}
\title{Print Parameter Selection Results in JPEX Deblurring}
\description{
  Display information about a JPEX parameter selection object.
}
\usage{
\method{print}{JPEX_Parameters}(x, type = "all", ...)
}
\arguments{
  \item{x}{A JPEX parameter selection object.}
  \item{type}{The type of information to display. The "cv_scores"
      option prints the cross-validation scores for each
      bandwidth. The "bandwidth" option prints the selected bandwidth.
      The "sigma" option prints the estimated noise level. The "all"
      option prints all the information.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Prints some information about a JPEX parameter selection object.
    In particular, this method prints the cross-validation scores,
    the selected bandwidth and estimated noise level.
}
\value{
    A display of parameter selection results in JPEX estimation.
}
\seealso{
    \code{\link{cv.jpex}}
}
\examples{
parSel <- cv.jpex(stopsign, c(2,3))
print(parSel, type = "all")
}
\author{Yicheng Kang}
\references{
Kang, Y. (2020) ``Consistent Blind Image Deblurring Using
Jump-Preserving Extrapolation'', \emph{Journal of Computational and
 Graphical Statistics}, \bold{29}(2), 372 -- 382,
 \doi{10.1080/10618600.2019.1665536}.
}
