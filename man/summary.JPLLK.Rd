\name{summary.JPLLK_Parameters}
\alias{summary.JPLLK_Parameters}
\title{Summarize Parameter Selection Results in JPLLK Estimation}
\description{
  Summarize and display some key information about a JPLLK parameter selection object.
}
\usage{
\method{summary}{JPLLK_Parameters}(object, ...)
}
\arguments{
  \item{object}{A JPLLK parameter selection object.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Summarize some information about a JPLLK parameter selection object.
    In particular, it displays the selected bandwidth and estimated noise
    level.
}
\value{
    A brief summary of parameter selection results in JPLLK estimation.
}
\seealso{
    \code{\link{JPLLK_surface}}, \code{\link{print.JPLLK_Parameters}},
    \code{\link{plot.JPLLK_Parameters}}
}
\examples{
fit <- JPLLK_surface(image=sar, bandwidth=c(3, 4))
summary(fit)
}
\author{Yicheng Kang}
\references{
Qiu, P. (2009) "Jump-Preserving Surface Reconstruction from Noisy Data",
\emph{Annals of the Institute of Statistical Mathematics}, \bold{61}(3),
715 -- 751, \doi{10.1007/s10463-007-0166-9}.
}
