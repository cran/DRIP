\name{plot.JPEX_Parameters}
\alias{plot.JPEX_Parameters}
\title{Plot Parameter Selection Results in JPEX Deblurring}
\description{
  Plot the results of a JPEX parameter selection procedure.
}
\usage{
\method{plot}{JPEX_Parameters}(x, ...)
}
\arguments{
  \item{x}{A JPEX parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Plot the cross-validation scores against the bandwidth values.
}
\value{
    A plot of cross-validation scores is produced. 
}
\seealso{
    \code{\link{cv.jpex}}, \code{\link{print.JPEX_Parameters}},
    \code{\link{summary.JPEX_Parameters}}
}
\examples{
parSel <- cv.jpex(stopsign, c(2,3))
plot(parSel)
}
\author{Yicheng Kang}
\references{
Kang, Y. (2020) ``Consistent Blind Image Deblurring Using
Jump-Preserving Extrapolation'', \emph{Journal of Computational and
 Graphical Statistics}, \bold{29}(2), 372 -- 382,
 \doi{10.1080/10618600.2019.1665536}.
}
