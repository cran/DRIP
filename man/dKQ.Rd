\name{dKQ}
\alias{dKQ}
\keyword{performance measure}
\title{edge detection, performance measure}
\description{ Compute the d_KQ distance between two sets of
  edge pixels. It can be used as a performance measure for
  (step/roof) edge detectors}
\usage{dKQ(edge1, edge2)}
\arguments{
  \item{edge1}{One set of pixels}
  \item{edge2}{The other set of pixels}
}
\value{Value of the \eqn{d_{KQ}}}
\details{
	The mathematical definition of \eqn{d_{KQ}} is as follows.
	\eqn{d_{KQ}(S_1, S_2) = \frac{0.5}{|S_1|}\sum_{p_1\in S_1} d_E(p_1, S_2) +
	\frac{0.5}{|S_2|}\sum_{p_2\in S_2} d_E(p_2, S_1)}, where \eqn{S_1} and \eqn{S_2}
	are two point sets, and \eqn{d_E} denotes the Euclidean distance.
}
\examples{
mat1 = matrix(c(1, rep(0, 3)), ncol = 2)
mat2 = matrix(c(rep(0, 3), 1), ncol = 2)
dKQ(mat1, mat2)
}

\references{
Kang, Y., and Qiu, P., "Jump Detection in Blurred Regression 
Surfaces," \emph{Technometrics}, \bold{56}, 2014, 539-550.
}