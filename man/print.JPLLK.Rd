\name{print.JPLLK_Parameters}
\alias{print.JPLLK_Parameters}
\title{Print Parameter Selection Results in JPLLK Estimation}
\description{
  Display information about a JPLLK parameter selection object.
}
\usage{
\method{print}{JPLLK_Parameters}(x, type =  "all", ...)
}
\arguments{
  \item{x}{A JPLLK parameter selection object.}
  \item{type}{The type of information to display. The "bandwidth" option
  prints the selected bandwidth. The "sigma" option prints the estimated
  noise level. The "all" option prints all the information.}
  \item{...}{Further arguments passed to or from other methods.}
}

\details{
    Prints some information about a JPLLK parameter selection object.
    In particular, this method prints the selected bandwidth and
    estimated noise level.
}
\value{
    A display of parameter selection results in JPLLK estimation.
}
\seealso{
    \code{\link{JPLLK_surface}}
}
\examples{
fit <- JPLLK_surface(image = sar, bandwidth = c(3, 4))
print(fit, type = "bandwidth")
}
\author{Yicheng Kang}
\references{
Qiu, P. (2009) "Jump-Preserving Surface Reconstruction from Noisy Data",
\emph{Annals of the Institute of Statistical Mathematics}, \bold{61}(3),
715 -- 751, \doi{10.1007/s10463-007-0166-9}.
}
