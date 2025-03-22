\name{posdefMatrix-class}
\title{Virtual Class \sQuote{posdefMatrix} of Positive Semidefinite Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{posdefMatrix-class}
%
\description{
  \code{posdefMatrix} is a virtual subclass of
  \code{\linkS4class{symmetricMatrix}} representing positive
  semidefinite matrices.  Positive semidefinite matrices tend to arise
  as the product of a matrix and its conjugate tranpose
  (\code{ct(x) \%*\% x} or \code{x \%*\% ct(x)}; see
  \code{\link{ct-methods}}).  In particular, covariance and
  correlation matrices are positive semidefinite.

  \code{posdefMatrix} was not defined until \pkg{Matrix} version 1.8-0.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{factors}}{inherited from virtual superclass
      \code{\linkS4class{symmetricMatrix}}.}
  }
}
\section{Subclasses}{
  The direct, nonvirtual subclasses of \code{posdefMatrix} defined in
  package \pkg{Matrix} are listed below with links to virtual
  superclasses defining their data type and storage format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Storage format\cr
    \code{dpoMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zpoMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dppMatrix} \tab \code{\linkS4class{dMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{zppMatrix} \tab \code{\linkS4class{zMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dpCMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zpCMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dpRMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zpRMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dpTMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zpTMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}
  }
  The subclasses of \code{zMatrix}, \code{CsparseMatrix},
  \code{RsparseMatrix}, and \code{TsparseMatrix} were not defined until
  \pkg{Matrix} version 1.8-0.

  Each of the above extends a nonvirtual class \code{ds.Matrix} or
  \code{zs.Matrix} representing Hermitian or symmetric matrices more
  generally; \code{dpoMatrix} and \code{dppMatrix} are extended by
  nonvirtual classes \code{corMatrix} and \code{copMatrix} representing
  correlation matrices.  For implementation details, see the help topic
  for the storage format.
}
\note{
  Tests for positive semidefiniteness in finite precision arithmetic are
  generally expensive and inexact.  Validity methods for nonvirtual
  subclasses of \code{posdefMatrix} defined in package \pkg{Matrix} test
  only that the diagonal entries of the represented matrix are
  nonnegative (a necessary, not sufficient condition).  Meanwhile,
  methods for \emph{coercion} to \code{posdefMatrix} defined in package
  \pkg{Matrix} test that the represented matrix has a Cholesky
  factorization (a sufficient, not necessary condition, being equivalent
  to positive \emph{definiteness}) by evaluating a suitable call to
  \code{\linkS4methods{Cholesky}}.  A more rigorous test for positive
  semidefiniteness, supporting matrices which are not positive definite,
  would use the Bunch-Kaufman or pivoted Cholesky factorization; see
  \code{\linkS4methods{BunchKaufman}} and
  \code{\linkS4methods{Cholesky}}.
}
\seealso{
  Related classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{triangularMatrix}},
  \code{\linkS4class{diagonalMatrix}}, and
  \code{\linkS4class{indexMatrix}} representing
  general (unstructured, possibly nonsquare), Hermitian or symmetric,
  triangular, diagonal, and index matrices.  Generic functions
  \code{\linkS4methods{crossprod}} and \code{\linkS4methods{tcrossprod}}
  for a common way to generate positive semidefinite matrices.  Generic
  function \code{\linkS4methods{Cholesky}} for (pivoted or unpivoted)
  Cholesky factorization of positive semidefinite matrices.  Generic
  function \code{\linkS4methods{chol}} for getting the upper triangular
  Cholesky factor.
}
\examples{
h6 <- Hilbert(6)
rcond(h6)
str(h6)
h6 * 27720 # is ``integer''
solve(h6)
str(hp6 <- pack(h6))

### Note that  as(*, "corMatrix")  *scales* the matrix
(ch6 <- as(h6, "corMatrix"))
stopifnot(all.equal(as(h6 * 27720, "dsyMatrix"), round(27720 * h6),
                    tolerance = 1e-14),
          all.equal(ch6@sd^(-2), 2*(1:6)-1,
                    tolerance = 1e-12))
chch <- Cholesky(ch6, perm = FALSE)
stopifnot(identical(chch, ch6@factors[["denseCholesky+-"]]),
          all(abs(crossprod(as(chch, "dtrMatrix")) - ch6) < 1e-10))
}
