\name{denseLU-class}
\title{Dense LU Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{denseLU-class}
%
\alias{coerce,denseLU,Matrix-method}
%
\description{
  \code{denseLU} is the class of dense, row-pivoted LU factorizations
  of real, square \emph{or non-square} matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{x}}{a numeric vector of length \code{prod(Dim)} storing
      the \code{L} and \code{U} factors together in a packed format.
      The details of the representation are specified by the LAPACK
      manual for its routine \code{DGETRF}.}
    \item{\code{perm}}{an integer vector of length \code{min(Dim)}
      specifying row interchanges as described in the LAPACK manual.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{LU}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{LU}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("denseLU", ...)}, but they are more typically obtained
  as the value of \code{\link{lu}(x)} for \code{x} inheriting from
  \code{\linkS4class{denseMatrix}} (often \code{\linkS4class{dgeMatrix}}).
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "denseLU", to = "Matrix")}:
      returns a \code{dtrMatrix} or \code{dtpMatrix}, useful for
      inspecting the internal representation of the factorization;
      see Note.}
    \item{expand}{\code{signature(x = "denseLU")}:
      see \code{\link{expand-methods}}.}
    \item{expand2}{\code{signature(x = "denseLU")}:
      see \code{\link{expand2-methods}}.}
    \item{solve}{\code{signature(a = "denseLU", b = "missing")}:
      see \code{\link{solve-methods}}.}
  }
}
\details{
  The decomposition is of the form
  \deqn{A = P L U}{A = P L U}
  where typically all matrices are of size \eqn{n\times n}{n by n}, and
  the matrix \eqn{P} is a permutation matrix, \eqn{L} is lower
  triangular and \eqn{U} is upper triangular (both of class
  \code{\linkS4class{dtrMatrix}}).

  Note that the \emph{dense} decomposition is also implemented for
  a \eqn{m\times n}{m by n} matrix \eqn{A}, when \eqn{m \ne n}{m != n}.

  If \eqn{m < n} (\dQuote{wide case}), \eqn{U} is \eqn{m\times n}{m by n},
  and hence not triangular.\cr
  If \eqn{m > n} (\dQuote{long case}), \eqn{L} is \eqn{m\times n}{m by n},
  and hence not triangular.
}
\seealso{
  Class \code{\linkS4class{sparseLU}} for sparse LU factorizations.
  
  Class \code{\linkS4class{dgeMatrix}}.

  Generic functions \code{\link{expand2}} and \code{\link{lu}}.
}
\examples{
set.seed(1)
mm <- Matrix(round(rnorm(9),2), nrow = 3)
mm
str(lum <- lu(mm))
elu <- expand(lum)
elu # three components: "L", "U", and "P", the permutation
elu$L \%*\% elu$U
(m2 <- with(elu, P \%*\% L \%*\% U)) # the same as 'mm'
stopifnot(all.equal(as(mm, "matrix"),
                    as(m2, "matrix")))

## See the underlying LAPACK representation
as(lum, "Matrix")
}
