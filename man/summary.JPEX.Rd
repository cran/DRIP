\name{summary.JPEX_Parameters}
\alias{summary.JPEX_Parameters}
\title{Summarize Parameter Selection Results in JPEX Deblurring}
\description{
  Summarize and display key information about a JPEX parameter selection object.
}
\usage{
\method{summary}{JPEX_Parameters}(object, ...)
}
\arguments{
  \item{object}{A JPEX parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Summarize some key information about a JPEX parameter selection object.
    In particular, it displays the chosen bandwidth and estimated noise
    level.
}
\value{
    A brief display of parameter selection results in JPEX estimation. 
}
\seealso{
    \code{\link{cv.jpex}}, \code{\link{print.JPEX_Parameters}},
    \code{\link{plot.JPEX_Parameters}}
}
\examples{
parSel <- cv.jpex(stopsign, c(2,3))
summary(parSel)
}
\author{Yicheng Kang}
\references{
Kang, Y. (2020) ``Consistent Blind Image Deblurring Using
Jump-Preserving Extrapolation'', \emph{Journal of Computational and
 Graphical Statistics}, \bold{29}(2), 372 -- 382,
 \doi{10.1080/10618600.2019.1665536}.
}
