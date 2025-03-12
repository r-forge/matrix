\name{denseMatrix-class}
\title{Virtual Class \sQuote{denseMatrix} of Dense Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{denseMatrix-class}
%
\description{
  \code{denseMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing dense format matrices.
  Formally, a dense format is a storage format in which the space
  required to represent an \eqn{m \times n}{m-by-n} matrix with \eqn{p}
  nonzero entries is not \eqn{O(\max\{m, n, p\})}{O(max(m, n, p))}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
  }
}
\section{Subclasses}{
  Package \pkg{Matrix} defines two virtual subclasses of
  \code{denseMatrix}: \code{\linkS4class{unpackedMatrix}} and
  \code{\linkS4class{packedMatrix}}.  These correspond to the
  conventional and packed storage formats defined by LAPACK:
  \url{https://netlib.org/lapack/lug/node121.html}.  See those help
  topics for implementation details and (recursively) further
  subclasses.
}
\seealso{
  \dQuote{Complementary} class \code{\linkS4class{sparseMatrix}}
  representing sparse (meaning \emph{not} dense) format matrices.
}
\examples{
showClass("denseMatrix")
}
