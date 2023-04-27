\name{expand-methods}
\title{Expand Matrix Factorizations}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{expand}
\alias{expand-methods}
\alias{expand2}
\alias{expand2-methods}
%
\alias{expand,CHMfactor-method}
\alias{expand,denseLU-method}
\alias{expand,sparseLU-method}
%
\alias{expand2,CHMfactor-method}
\alias{expand2,Cholesky-method}
\alias{expand2,Schur-method}
\alias{expand2,denseLU-method}
\alias{expand2,pCholesky-method}
\alias{expand2,sparseLU-method}
%
\description{
  \code{expand} and \code{expand2} construct matrix factors from
  objects specifying matrix factorizations.  Such objects typically
  do not store the factors explicitly, employing instead a compact
  representation in order to save memory.
}
\usage{
expand2(x, \dots)
expand (x, \dots)
}
\arguments{
  \item{x}{a matrix factorization, typically inheriting from
    virtual class \code{\linkS4class{MatrixFactorization}}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  A list of matrix factors, typically with names using conventional
  notation, as in \code{list(L=, U=)}.
}
\details{
  Methods for \code{expand} are retained only for backwards
  compatibility with \pkg{Matrix} \code{<= 1.5-4}.

  New code should use \code{expand2} whose methods are more consistent,
  all obeying the rule that the product of the returned matrix factors
  should reproduce (within some tolerance) the factorized matrix,
  \emph{including} its \code{dimnames}.  Hence if \code{x} is a matrix
  and \code{y} is its factorization, then
  \preformatted{all.equal(as(x, "matrix"), as(Reduce(`\%*\%`, expand2(y)), "matrix"))}
  should in most cases return \code{TRUE}.
}
\section{Methods}{
  \describe{
    \item{expand2}{\code{signature(x = "CHMfactor")}:
      lists the factors in one of the factorizations
      \eqn{A = P' L L' P = P' L_{1} D L_{1}' P}
      of \eqn{n}-by-\eqn{n} \eqn{A}, where
      \eqn{P} is a permutation matrix
      (of class \code{\linkS4class{pMatrix}}),
      \eqn{D} is a diagonal matrix
      (of class \code{\linkS4class{ddiMatrix}}),
      \eqn{L} and \eqn{L_{1}} are nonunit and unit lower triangular matrices
      (of class \code{\linkS4class{dtCMatrix}})
      related by \eqn{L = L_{1} sqrt(D)}.
      Which factorization depends on optional argument \code{LDL},
      whose default value is \code{TRUE}.  It does \emph{not} depend
      on \code{isLDL(x)}.
      The names of the result are \code{c("P.", "L", "L.", "P")} or
      \code{c("P.", "L1", "D", "L1.", "P")}.}
    \item{expand2}{\code{signature(x = "Cholesky")}:
      lists the factors in one of the factorizations
      \eqn{A = L L' = L_{1} D L_{1}'}
      of \eqn{n}-by-\eqn{n} \eqn{A}, where
      \eqn{D} is a diagonal matrix
      (of class \code{\linkS4class{ddiMatrix}}) and 
      \eqn{L} and \eqn{L_{1}} are nonunit and unit lower triangular matrices
      (of class \code{\linkS4class{dtrMatrix}})
      related by \eqn{L = L_{1} sqrt(D)}.
      Which factorization depends on optional argument \code{LDL},
      whose default value is \code{TRUE}.
      The names of the result are \code{c("L", "L.")} or
      \code{c("L1", "D", "L1.")}.}
    \item{expand2}{\code{signature(x = "pCholesky")}:
      similar to \code{signature(x = "Cholesky")},
      except that \eqn{L} and \eqn{L_{1}} are constructed as
      \code{\linkS4class{dtpMatrix}}.}
    \item{expand2}{\code{signature(x = "denseLU")}:
      lists the factors in the factorization \eqn{A = P L U}
      of \eqn{m}-by-\eqn{n} \eqn{A}, where
      \eqn{P} is a permutation matrix
      (of class \code{\linkS4class{pMatrix}}),
      \eqn{L} is a unit lower trapezoidal matrix
      (of class \code{\linkS4class{dtrMatrix}} if \eqn{m <= n},
      \code{\linkS4class{dgeMatrix}} otherwise), and
      \eqn{L} is an upper trapezoidal matrix
      (of class \code{\linkS4class{dtrMatrix}} if \eqn{m >= n},
      \code{\linkS4class{dgeMatrix}} otherwise).
      The names of the result are \code{c("P", "L", "U")}.}
    \item{expand2}{\code{signature(x = "sparseLU")}:
      lists the factors in the factorization \eqn{A = P' L U Q}
      of \eqn{n}-by-\eqn{n} \eqn{A}, where
      \eqn{P} and \eqn{Q} are permutation matrices
      (of class \code{\linkS4class{pMatrix}}),
      \eqn{L} is a unit lower triangular matrix
      (of class \code{\linkS4class{dtCMatrix}}), and
      \eqn{U} is an upper triangular matrix
      (of class \code{\linkS4class{dtCMatrix}}).
      The names of the result are \code{c("P.", "L", "U", "Q")}.}
    \item{expand2}{\code{signature(x = "Schur")}:
      lists the factors in the factorization \eqn{A = Q T Q'}
      of \eqn{n}-by-\eqn{n} \eqn{A}, where
      \eqn{Q} is a orthogonal matrix
      (of virtual class \code{\linkS4class{dMatrix}})
      and \eqn{T} is a block upper triangular matrix
      (of virtual class \code{\linkS4class{dMatrix}})
      with 1-by-1 and 2-by-2 diagonal blocks.
      The names of the result are \code{c("Q", "T", "Q.")}.}
    \item{expand}{\code{signature(x = "CHMfactor")}:
      as \code{expand2}, but listing only the factors
      \eqn{P} and \eqn{L} in \eqn{A = P' L L' P}.
      The result has names \code{c("P", "L")},
      and all of the components have empty \code{Dimnames} slots.}
    \item{expand}{\code{signature(x = "denseLU")}:
      as \code{expand2}, but listing factor \eqn{P}
      in \eqn{A = P L U} after \eqn{L} and \eqn{U}.
      The result has names \code{c("L", "U", "P")},
      and all of the components have empty \code{Dimnames} slots.}
    \item{expand}{\code{signature(x = "sparseLU")}:
      as \code{expand2}, but listing \eqn{P} instead of
      factor \eqn{P'} in \eqn{A = P' L U Q}.
      The result has names \code{c("P", "L", "U", "Q")}.
      Components \code{L} and \code{U} store the row and column names,
      of the factorized matrix in their \code{Dimnames} slots.
      Components \code{P} and \code{Q} have empty \code{Dimnames} slots.}
  }
}
\seealso{
  The virtual class \code{\linkS4class{compMatrix}}
  of \emph{factorizable} matrices.

  The virtual class \code{\linkS4class{MatrixFactorization}}
  of matrix factorizations.
  
  Generic functions \code{\link{chol}}, \code{\link{Cholesky}},
  \code{\link{BunchKaufman}}, \code{\link{Schur}}, \code{\link{lu}},
  and \code{\link{qr}} for \emph{computing} factorizations.
}
\examples{
(x <- Matrix(round(rnorm(9),2), 3, 3))
(ex <- expand(lux <- lu(x)))
}
