\name{expm}
\alias{expm}
\alias{expm,dgeMatrix-method}
\title{Matrix exponential}
\description{
  Return the exponential of a matrix.
  
}
\usage{
expm(x)
}
\arguments{
  \item{x}{An object that inherits from the \code{Matrix} class}
}
\details{
  The exponential of a matrix is defined as the infinite Taylor
  series \code{expm(A) = I + A + A^2/2! + A^3/3! + ...} (although this is
  definitely not the way to compute it).  The method for the
  \code{dgeMatrix} class uses Ward's diagonal Pade' approximation with
  three step preconditioning.  It is a simple translation of the
  implementation of the corresponding Octave function contributed
  to the Octave project by A. Scottedward Hodel
  \email{A.S.Hodel@Eng.Auburn.EDU}.
}
\value{
  The matrix exponential of \code{x}.
}
\references{}
%\author{}
%\note{}
%\seealso{}
\examples{
(mm <- Matrix(c(-49, -64, 24, 31), nc = 2))
expm(mm)
(mm <- Matrix(cbind(0,rbind(6*diag(3),0)), nc = 4))
expm(mm)
(mm <- Matrix(c(1,0,1,1), nc = 2))
expm(mm)
}
\keyword{algebra}
\keyword{math}
