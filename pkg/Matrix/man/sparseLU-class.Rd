\name{sparseLU-class}
\title{Sparse LU Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{sparseLU-class}
%
\alias{determinant,sparseLU,logical-method}
%
\description{
  \code{sparseLU} is the class of sparse, row- and column-pivoted
  LU factorizations of \eqn{n \times n}{n-by-n} real matrices \eqn{A},
  having the general form
  \deqn{P_{1} A P_{2} = L U}{P1 * A * P2 = L * U}
  or (equivalently)
  \deqn{A = P_{1}' L U P_{2}'}{A = P1' * L * U * P2'}
  where
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices,
  \eqn{L} is a unit lower triangular matrix, and
  \eqn{U} is an upper triangular matrix.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{L}}{an object of class \code{\linkS4class{dtCMatrix}},
      the unit lower triangular \eqn{L} factor.}
    \item{\code{U}}{an object of class \code{\linkS4class{dtCMatrix}},
      the upper triangular \eqn{U} factor.}
    \item{\code{p}, \code{q}}{0-based integer vectors of length
      \code{Dim[1]},
      specifying the permutations applied to the rows and columns of
      the factorized matrix.  \code{q} of length 0 is valid and
      equivalent to the identity permutation, implying no column pivoting.
      Using \R{} syntax, the matrix \eqn{P_{1} A P_{2}}{P1 * A * P2}
      is precisely \code{A[p+1, q+1]}
      (\code{A[p+1, ]} when \code{q} has length 0).}
  }
}
\section{Extends}{
  Class \code{\linkS4class{LU}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{LU}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("sparseLU", ...)}, but they are more typically obtained
  as the value of \code{\link{lu}(x)} for \code{x} inheriting from
  \code{\linkS4class{sparseMatrix}} (often \code{\linkS4class{dgCMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "sparseLU", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand}}{\code{signature(x = "sparseLU")}:
      see \code{\link{expand-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "sparseLU")}:
      see \code{\link{expand2-methods}}.}
  }
}
\seealso{
  Class \code{\linkS4class{denseLU}} for dense LU factorizations.
  
  Class \code{\linkS4class{dgCMatrix}}.

  Generic functions \code{\link{expand2}} and \code{\link{lu}}.
}
\references{
  Davis, T. A. (2006).
  \emph{Direct Methods for Sparse Linear Systems},
  Society for Industrial and Applied Mathematics.
}
\examples{
## Extending the one in   examples(lu), calling the matrix  A,
## and confirming the factorization identities :
A <- as(readMM(system.file("external/pores_1.mtx",
                            package = "Matrix")),
         "CsparseMatrix")
## with dimnames(.) - to see that they are propagated :
dimnames(A) <- list(paste0("r", seq_len(nrow(A))),
                    paste0("C", seq_len(ncol(A))))
str(luA <- lu(A)) # p is a 0-based permutation of the rows
                  # q is a 0-based permutation of the columns
xA <- expand(luA)
stopifnot(identical(unname(xA$L), luA@L),
          identical(unname(xA$U), luA@U),
          identical(xA$P, as(luA@p + 1L, "pMatrix")),
          identical(xA$Q, as(luA@q + 1L, "pMatrix")))

P.LUQ <- with(xA, t(P) \%*\% L \%*\% U \%*\% Q)
stopifnot(all.equal(unname(A), unname(P.LUQ), tolerance = 1e-12))

## permute rows and columns of original matrix
pA <- A[luA@p + 1L, luA@q + 1L]
PAQ. <- with(xA, P \%*\% A \%*\% t(Q))
stopifnot(all.equal(unname(pA), unname(PAQ.), tolerance = 1e-12))

pLU <- drop0(xA$L \%*\% xA$U) # L \%*\% U -- dropping extra zeros
stopifnot(all.equal(pA, pLU, tolerance = 1e-12)) # (incl. permuted row- and column-names)
}
