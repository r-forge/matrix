\name{sparseQR-class}
\title{Sparse QR Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
\keyword{utilities}
%
\alias{sparseQR-class}
%
\alias{determinant,sparseQR,logical-method}
\alias{qr.Q,sparseQR-method}
\alias{qr.R,sparseQR-method}
\alias{qr.X,sparseQR-method}
\alias{qr.coef,sparseQR,Matrix-method}
\alias{qr.coef,sparseQR,dgeMatrix-method}
\alias{qr.coef,sparseQR,matrix-method}
\alias{qr.coef,sparseQR,numLike-method}
\alias{qr.fitted,sparseQR,Matrix-method}
\alias{qr.fitted,sparseQR,dgeMatrix-method}
\alias{qr.fitted,sparseQR,matrix-method}
\alias{qr.fitted,sparseQR,numLike-method}
\alias{qr.qty,sparseQR,Matrix-method}
\alias{qr.qty,sparseQR,dgeMatrix-method}
\alias{qr.qty,sparseQR,matrix-method}
\alias{qr.qty,sparseQR,numLike-method}
\alias{qr.qy,sparseQR,Matrix-method}
\alias{qr.qy,sparseQR,dgeMatrix-method}
\alias{qr.qy,sparseQR,matrix-method}
\alias{qr.qy,sparseQR,numLike-method}
\alias{qr.resid,sparseQR,Matrix-method}
\alias{qr.resid,sparseQR,dgeMatrix-method}
\alias{qr.resid,sparseQR,matrix-method}
\alias{qr.resid,sparseQR,numLike-method}
%
\alias{qrR}
%
\description{
  \code{sparseQR} is the class of sparse, row- and column-pivoted
  QR factorizations of \eqn{m \times n}{m-by-n} (\eqn{m \ge n}{m >= n})
  real matrices, having the general form
  \deqn{P_{1} A P_{2} = Q R = \begin{bmatrix} Q_{1} & Q_{2} \end{bmatrix} \begin{bmatrix} R_{1} \\ 0 \end{bmatrix} = Q_{1} R_{1}}{P1 * A * P2 = Q * R = [Q1, Q2] * [R1; 0] = Q1 * R1}
  or (equivalently)
  \deqn{A = P_{1}' Q R P_{2}' = P_{1}' \begin{bmatrix} Q_{1} & Q_{2} \end{bmatrix} \begin{bmatrix} R_{1} \\ 0 \end{bmatrix} P_{2}' = P_{1}' Q_{1} R_{1} P_{2}'}{A = P1' * Q * R * P2' = P1' * [Q1, Q2] * [R1; 0] * P2' = P1' * Q1 * R1 * P2'}
  where
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices,
  \eqn{Q = \prod_{j = 1}^{n} Q_{j}}{Q = prod(Qj : j = 1,...,n)}
  is an \eqn{m \times m}{m-by-m} orthogonal matrix
  (\eqn{Q_{1}}{Q1} contains the first \eqn{n} column vectors)
  equal to the product of \eqn{n} Householder matrices \eqn{Q_{j}}{Qj}, and
  \eqn{R} is an \eqn{m \times n}{m-by-n} upper trapezoidal matrix
  (\eqn{R_{1}}{R1} contains the first \eqn{n} row vectors and is
  upper \emph{triangular}) with non-negative diagonal elements.
}
\usage{
qrR(qr, complete = FALSE, backPermute = TRUE, row.names = TRUE)
}
\arguments{
  \item{qr}{an object of class \code{\linkS4class{sparseQR}},
    almost always the result of a call to generic function \code{qr}
    with sparse \code{x}.}
  \item{complete}{a logical indicating if \eqn{R} should be returned
    instead of \eqn{R_{1}}{R1}.  Note that \eqn{R} and \eqn{R_{1}}{R1}
    coincide in the \eqn{m = n} case.}
  \item{backPermute}{a logical indicating if \eqn{R_{1} P_{2}'}{R1 * P2'} 
    (\eqn{R P_{2}'}{R * P2'} if \code{complete = TRUE}) should be returned 
    instead of \eqn{R_{1}}{R1} (\eqn{R} if \code{complete = TRUE}).}
  \item{row.names}{a logical indicating if \code{dimnames(qr)[1]}
    should be propagated unpermuted to the result.
    If \code{complete = FALSE}, then only the first \eqn{n} names are kept.}
}
\value{
  \code{qrR} returns
  a \code{\linkS4class{dtCMatrix}} representing \eqn{R_{1}}{R1}
  (if \code{complete = FALSE}) or
  a \code{\linkS4class{dgCMatrix}} representing \eqn{R}
  (if \code{complete = TRUE}).
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{beta}}{a numeric vector of length \code{Dim[2]},
      used to construct Householder matrices; see \code{V} below.}
    \item{\code{V}}{an object of class \code{\linkS4class{dgCMatrix}}
      with \code{Dim[2]} columns.  The number of rows \code{nrow(V)}
      is at least \code{Dim[1]} and at most \code{Dim[1]+Dim[2]}.
      \code{V} is lower trapezoidal, and its column vectors generate the
      Householder matrices \eqn{Q_{j}}{Qj} that compose the orthogonal
      \eqn{Q} factor.  Specifically, \eqn{Q_{j}}{Qj} is constructed as
      \code{diag(Dim[1]) - beta[j] * tcrossprod(V[, j])}.}
    \item{\code{R}}{an object of class \code{\linkS4class{dgCMatrix}}
      with \code{nrow(V)} rows and \code{Dim[2]} columns.
      \code{R} is the upper trapezoidal \eqn{R} factor with
      non-negative diagonal elements.}
    \item{\code{p}, \code{q}}{0-based integer vectors of length
      \code{nrow(V)} and \code{Dim[2]}, respectively,
      specifying the permutations applied to the rows and columns of
      the factorized matrix.  \code{q} of length 0 is valid and
      equivalent to the identity permutation, implying no column pivoting.
      Using \R{} syntax, the matrix \eqn{P_{1} A P_{2}}{P1 * A * P2}
      is precisely \code{A[p+1, q+1]}
      (\code{A[p+1, ]} when \code{q} has length 0).}
  }
}
\section{Extends}{
  Class \code{\linkS4class{QR}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{QR}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("sparseQR", ...)}, but they are more typically obtained
  as the value of \code{\link{qr}(x)} for \code{x} inheriting from
  \code{\linkS4class{sparseMatrix}} (often \code{\linkS4class{dgCMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "sparseQR", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand2}}{\code{signature(x = "sparseQR")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{qr.Q}}{\code{signature(qr = "sparseQR")}:
      returns as a \code{\linkS4class{dgeMatrix}} either
      \eqn{P_{1}' Q}{P1' * Q} or \eqn{P_{1}' Q_{1}}{P1' * Q1},
      depending on optional argument \code{complete}.  The default
      is \code{FALSE}, indicating \eqn{P_{1}' Q_{1}}{P1' * Q1}.}
    \item{\code{qr.R}}{\code{signature(qr = "sparseQR")}:
      returns either \eqn{R} as a \code{\linkS4class{dgCMatrix}}
      or \eqn{R_{1}}{R1} as a \code{\linkS4class{dtCMatrix}},
      depending on optional argument \code{complete}.  The default
      is \code{FALSE}, indicating \eqn{R_{1}}{R1}.}
    \item{\code{qr.X}}{\code{signature(qr = "sparseQR")}:
      returns \eqn{A} as a \code{\linkS4class{dgeMatrix}},
      by default.  If \eqn{m > n} and optional argument
      \code{ncol} is greater than \eqn{n}, then the result
      is augmented with \eqn{P_{1}' Q J}{P1 * Q * J}, where
      \eqn{J} is composed of columns \eqn{(n+1)} through
      \code{ncol} of the \eqn{m \times m}{m-by-m} identity matrix.}
    \item{\code{qr.coef}}{\code{signature(qr = "sparseQR", y = .)}:
      returns as a \code{\linkS4class{dgeMatrix}} or vector
      the result of multiplying \code{y} on the left by
      \eqn{P_{2} R_{1}^{-1} Q_{1}' P_{1}}{P2 * R1^{-1} * Q1' * P1}.}
    \item{\code{qr.fitted}}{\code{signature(qr = "sparseQR", y = .)}:
      returns as a \code{\linkS4class{dgeMatrix}} or vector
      the result of multiplying \code{y} on the left by
      \eqn{P_{1}' Q_{1} Q_{1}' P_{1}}{P1' * Q1 * Q1' * P1}.}
    \item{\code{qr.resid}}{\code{signature(qr = "sparseQR", y = .)}:
      returns as a \code{\linkS4class{dgeMatrix}} or vector
      the result of multiplying \code{y} on the left by
      \eqn{P_{1}' Q_{2} Q_{2}' P_{1}}{P1' * Q2 * Q2' * P1}.}
    \item{\code{qr.qty}}{\code{signature(qr = "sparseQR", y = .)}:
      returns as a \code{\linkS4class{dgeMatrix}} or vector
      the result of multiplying \code{y} on the left by
      \eqn{Q' P_{1}}{Q' * P1}.}
    \item{\code{qr.qy}}{\code{signature(qr = "sparseQR", y = .)}:
      returns as a \code{\linkS4class{dgeMatrix}} or vector
      the result of multiplying \code{y} on the left by
      \eqn{P_{1}' Q}{P1' * Q}.}
    \item{\code{solve}}{\code{signature(a = "sparseQR", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\details{
  The method for \code{qr.Q} does not return \eqn{Q} but rather the
  (also orthogonal) product \eqn{P_{1}' Q}{P1' * Q}.  This behaviour
  is algebraically consistent with the \pkg{base} implementation
  (see \code{\link[base]{qr}}), which can be seen by noting that
  \code{qr.default} in \pkg{base} does not pivot rows, constraining
  \eqn{P_{1}}{P1} to be an identity matrix.  It follows that
  \code{qr.Q(qr.default(x))} also returns \eqn{P_{1}' Q}{P1' * Q}.
  
  Similarly, the methods for \code{qr.qy} and \code{qr.qty} multiply
  on the left by \eqn{P_{1}' Q}{P1' * Q} and \eqn{Q' P_{1}}{Q' * P1}
  rather than \eqn{Q} and \eqn{Q'}.
  
  It is wrong to expect the values of \code{qr.Q} (or \code{qr.R},
  \code{qr.qy}, \code{qr.qty}) computed from \dQuote{equivalent}
  sparse and dense factorizations
  (say, \code{qr(x)} and \code{qr(as(x, "matrix"))} for \code{x}
  of class \code{\linkS4class{dgCMatrix}}) to compare equal.
  The underlying factorization algorithms are quite different,
  notably as they employ different pivoting strategies,
  and in general the factorization is not unique even for fixed
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2}.
  
  On the other hand, the values of \code{qr.X}, \code{qr.coef},
  \code{qr.fitted}, and \code{qr.resid} are well-defined, and
  in those cases the sparse and dense computations \emph{should}
  compare equal (within some tolerance).
}
\seealso{
  Class \code{\linkS4class{dgCMatrix}}.
  
  Generic function \code{\link[base]{qr}} from \pkg{base},
  whose default method \code{qr.default} \dQuote{defines}
  the S3 class \code{qr} of dense QR factorizations.

  \code{\link{qr-methods}} for methods defined in \pkg{Matrix}.
  
  The many auxiliary functions for QR factorizations:
  \code{\link{qr.Q}}, \code{\link{qr.R}}, \code{\link{qr.X}},
  \code{\link{qr.coef}}, \code{\link{qr.fitted}}, \code{\link{qr.resid}},
  \code{\link{qr.qty}}, \code{\link{qr.qy}}, and \code{\link{qr.solve}}.
}
\references{
  Davis, T. A. (2006).
  \emph{Direct Methods for Sparse Linear Systems},
  Society for Industrial and Applied Mathematics.
}
\examples{
showClass("sparseQR")

data(KNex, package = "Matrix")
A <- KNex[["mm"]]

## With dimnames, to see that they are propagated :
d <- dim(A)
dimnames(A) <- dn <- list(paste0("r", seq_len(d[1L])),
                          paste0("c", seq_len(d[2L])))

(qr.A <- qr(A))
e.qr.A <- expand2(qr.A)
E.qr.A <- expand2(qr.A, complete = TRUE)

d
t(sapply(e.qr.A, dim))
t(sapply(E.qr.A, dim))

## Horribly inefficient, but instructive :
slowQ <- function(V, beta) {
    d <- dim(V)
    Q <- diag(d[1L])
    if(d[2L] > 0L) {
        for(j in d[2L]:1L) {
            cat(j, "\n", sep = "")
            Q <- Q - (beta[j] * tcrossprod(V[, j])) \%*\% Q
        }
    }
    Q
}

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

## A ~ P1' Q R P2' ~ P1' Q1 R1 P2' in floating point
stopifnot(exprs = {
    identical(names(e.qr.A), c("P1.", "Q1", "R1", "P2."))
    identical(names(E.qr.A), c("P1.", "Q" , "R" , "P2."))
    identical(e.qr.A[["P1."]],
              new("pMatrix", Dim = d[c(1L, 1L)], Dimnames = c(dn[1L], list(NULL)),
                  margin = 1L, perm = invertPerm(qr.A@p, 0L, 1L)))
    identical(e.qr.A[["P2."]],
              new("pMatrix", Dim = d[c(2L, 2L)], Dimnames = c(list(NULL), dn[2L]),
                  margin = 2L, perm = invertPerm(qr.A@q, 0L, 1L)))
    identical(e.qr.A[["R1"]], triu(E.qr.A[["R"]][seq_len(d[2L]), ]))
    identical(e.qr.A[["Q1"]],      E.qr.A[["Q"]][, seq_len(d[2L])] )
    identical(E.qr.A[["R"]], qr.A@R)
 ## ae1(E.qr.A[["Q"]], slowQ(qr.A@V, qr.A@beta))
    ae1(crossprod(E.qr.A[["Q"]]), diag(d[1L]))
    ae1(A, with(e.qr.A, P1. \%*\% Q1 \%*\% R1 \%*\% P2.))
    ae1(A, with(E.qr.A, P1. \%*\% Q  \%*\% R  \%*\% P2.))
    ae2(A[qr.A@p + 1L, qr.A@q + 1L], with(e.qr.A, Q1 \%*\% R1))
    ae2(A[qr.A@p + 1L, qr.A@q + 1L], with(E.qr.A, Q  \%*\% R ))
})

## More identities
b <- KNex[["y"]]
B <- cbind(b, sin(b), deparse.level = 0)
stopifnot(exprs = {
    ae1(qr.X     (qr.A   ), A)
    ae2(qr.Q     (qr.A   ), with(e.qr.A, P1. \%*\% Q1))
    ae2(qr.R     (qr.A   ), with(e.qr.A, R1))
    ae2(qr.coef  (qr.A, B), with(e.qr.A, solve(R1 \%*\% P2., t(P1. \%*\% Q1)) \%*\% B))
    ae2(qr.fitted(qr.A, B), with(e.qr.A,     tcrossprod(P1. \%*\% Q1) \%*\% B))
    ae2(qr.resid (qr.A, B), with(e.qr.A, B - tcrossprod(P1. \%*\% Q1) \%*\% B))
    ae2(qr.qy    (qr.A, B), with(E.qr.A,   P1. \%*\% Q  \%*\% B))
    ae2(qr.qty   (qr.A, B), with(E.qr.A, t(P1. \%*\% Q) \%*\% B))
    ae2(qr.qy(qr.A, qr.qty(qr.A, B)), B)
})

## Sparse and dense computations should agree here
qr.Am <- qr(as(A, "matrix")) # <=> qr.default(A)
stopifnot(exprs = {
    ae2(qr.X     (qr.A   ), qr.X     (qr.Am   ))
    ae2(qr.coef  (qr.A, B), qr.coef  (qr.Am, B))
    ae2(qr.fitted(qr.A, B), qr.fitted(qr.Am, B))
    ae2(qr.resid (qr.A, B), qr.resid (qr.Am, B))
})
}
