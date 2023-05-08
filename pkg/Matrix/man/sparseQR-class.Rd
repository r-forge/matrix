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
\alias{qr.Q,sparseQR-method}
\alias{qr.R,sparseQR-method}
\alias{qr.coef,sparseQR,Matrix-method}
\alias{qr.coef,sparseQR,ddenseMatrix-method}
\alias{qr.coef,sparseQR,matrix-method}
\alias{qr.coef,sparseQR,numeric-method}
\alias{qr.fitted,sparseQR,Matrix-method}
\alias{qr.fitted,sparseQR,ddenseMatrix-method}
\alias{qr.fitted,sparseQR,matrix-method}
\alias{qr.fitted,sparseQR,numeric-method}
\alias{qr.qty,sparseQR,Matrix-method}
\alias{qr.qty,sparseQR,ddenseMatrix-method}
\alias{qr.qty,sparseQR,matrix-method}
\alias{qr.qty,sparseQR,numeric-method}
\alias{qr.qy,sparseQR,Matrix-method}
\alias{qr.qy,sparseQR,ddenseMatrix-method}
\alias{qr.qy,sparseQR,matrix-method}
\alias{qr.qy,sparseQR,numeric-method}
\alias{qr.resid,sparseQR,Matrix-method}
\alias{qr.resid,sparseQR,ddenseMatrix-method}
\alias{qr.resid,sparseQR,matrix-method}
\alias{qr.resid,sparseQR,numeric-method}
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
  \eqn{Q} is an \eqn{m \times m}{m-by-m} orthogonal matrix
  (\eqn{Q_{1}}{Q1} contains the first \eqn{n} column vectors), and
  \eqn{R} is an \eqn{m \times n}{m-by-n} upper trapezoidal matrix
  (\eqn{R_{1}}{R1} contains the first \eqn{n} row vectors and is
  upper \emph{triangular}).
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
  \item{row.names}{a logical indicating if \code{rownames(qr)}
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
    \item{\code{beta}}{a numeric vector of length \code{Dim[2]}.
      \code{beta} contains the normalizing constants used to compute
      the Householder matrices that compose the \eqn{Q} factor.}
    \item{\code{V}}{an object of class \code{\linkS4class{dgCMatrix}}
      with at least \code{Dim[1]} rows and exactly \code{Dim[2]} columns.  
      The column vectors of \code{V} generate the Householder matrices
      that compose the \eqn{Q} factor.}
    \item{\code{R}}{an object of class \code{\linkS4class{dgCMatrix}}
      with \code{nrow(V)} rows and \code{Dim[2]} columns.
      \code{R} is the upper trapezoidal \eqn{R} factor.}
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
      is \code{FALSE}, indicating \eqn{R_{1}}{R1}.
      This method warns about possible differences in the result of
      \code{qr.R(qr(x))} depending on whether \code{x} is sparse or
      dense.  The warning can be suppressed with
      \code{\link{options}(Matrix.quiet.qr.R = TRUE)}
      or
      \code{\link{options}(Matrix.quiet = TRUE)},
      the latter having side effects not limited to this method.}
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
  The orthogonal factor \eqn{Q} is defined implicitly by slots
  \code{V} and \code{beta}.  The columns of \code{V} and elements
  of \code{beta} generate the \eqn{n} Householder matrices whose
  product is \eqn{Q}.

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
  
  It is wrong to expect equality of values of \code{qr.Q}
  (or \code{qr.R}, \code{qr.qy}, \code{qr.qty}) computed from
  \dQuote{equivalent} sparse and dense factorizations (say,
  \code{qr(x)} and \code{qr(as(x, "matrix"))} for \code{x}
  of class \code{\linkS4class{dgCMatrix}}).  The underlying
  factorization algorithms are quite different (notably as
  they employ different criteria to select \eqn{P_{2}}{P2}),
  and in general the factorization is not unique even for
  fixed \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2}.
  
  On the other hand, the values of \code{qr.coef}, \code{qr.fitted}
  (with \code{k = ncol(qr@R)}), and \code{qr.resid} are well-defined,
  so the sparse and dense computations should agree.
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
data(KNex)
mm <- KNex $ mm
 y <- KNex $  y
 y. <- as(y, "CsparseMatrix")
str(qrm <- qr(mm))
 qc  <- qr.coef  (qrm, y); qc. <- qr.coef  (qrm, y.) # 2nd failed in Matrix <= 1.1-0
 qf  <- qr.fitted(qrm, y); qf. <- qr.fitted(qrm, y.)
 qs  <- qr.resid (qrm, y); qs. <- qr.resid (qrm, y.)
stopifnot(all.equal(qc, as.numeric(qc.),  tolerance=1e-12),
          all.equal(qf, as.numeric(qf.),  tolerance=1e-12),
          all.equal(qs, as.numeric(qs.),  tolerance=1e-12),
          all.equal(qf+qs, y, tolerance=1e-12))
%% FIXME: add more similarly to ./lu.Rd, see also ./qr-methods.Rd
}
