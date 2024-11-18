\name{plot.JPLLK_Parameters}
\alias{plot.JPLLK_Parameters}
\title{Plot the Residual Image from a JPLLK Estimation}
\description{
  Plot the residual image from a JPLLK parameter selection object where
  the bandwidth is chosen by cross validation.
}
\usage{
\method{plot}{JPLLK_Parameters}(x, ...)
}
\arguments{
  \item{x}{A JPLLK parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Plot the residual image from a JPLLK parameter selection object.
    The bandwidth used in producing the residual image is chosen by
    cross validation. 
}
\value{
    A residual image is produced. 
}
\seealso{
    \code{\link{JPLLK_surface}}, \code{\link{print.JPLLK_Parameters}},
    \code{\link{summary.JPLLK_Parameters}}
}
\examples{
fit <- JPLLK_surface(image=sar, bandwidth=c(3, 4))
plot(fit)
}
\author{Yicheng Kang}
\references{
Qiu, P. (2009) "Jump-Preserving Surface Reconstruction from Noisy Data",
\emph{Annals of the Institute of Statistical Mathematics}, \bold{61}(3),
715 -- 751, \doi{10.1007/s10463-007-0166-9}.
}
