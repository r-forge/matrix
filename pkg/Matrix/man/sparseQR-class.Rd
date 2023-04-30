\name{sparseQR-class}
\title{Sparse QR Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
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
    \item{\code{p}, \code{q}}{integer vectors of length \code{nrow(V)}
      and \code{Dim[2]} specifying row and column permutations,
      respectively.  \code{q} may have length 0, implying no column
      permutation.}
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
      returns the orthogonal matrix \eqn{P_{1}' Q}{P1' * Q}.}
    \item{\code{qr.R}}{\code{signature(qr = "sparseQR")}:
      returns the upper trapezoidal matrix \eqn{R} or the result
      \eqn{R P_{2}'}{R * P2'} of permuting its columns, the latter
      if \code{backPermute = TRUE} is passed in the call.
      This method warns about possible differences between its
      result and that of \code{base::qr.R(base::qr.default(.))}.
      The warning can be suppressed with
      \code{\link{options}(Matrix.quiet.qr.R = TRUE)}
      or
      \code{\link{options}(Matrix.quiet = TRUE)},
      the latter having side effects not limited to this method.}
    \item{\code{qr.coef}}{\code{signature(qr = "sparseQR", y = .)}:
      returns the result of multiplying \code{y} on the left
      by \eqn{P_{2} R_{1}^{-1} Q_{1}' P_{1}}{P2 * R1^{-1} * Q1' * P1}.}
    \item{\code{qr.fitted}}{\code{signature(qr = "sparseQR", y = .)}:
      returns the result of multiplying \code{y} on the left
      by \eqn{P_{1}' Q_{1} Q_{1}' P_{1}}{P1' * Q1 * Q1' * P1}.}
    \item{\code{qr.resid}}{\code{signature(qr = "sparseQR", y = .)}:
      returns the result of multiplying \code{y} on the left
      by \eqn{P_{1}' Q_{2} Q_{2}' P_{1}}{P1' * Q2 * Q2' * P1}.}
    \item{\code{qr.qty}}{\code{signature(qr = "sparseQR", y = .)}:
      returns the result of multiplying \code{y} on the left
      by \eqn{Q' P_{1}}{Q' * P1}.}
    \item{\code{qr.qy}}{\code{signature(qr = "sparseQR", y = .)}:
      returns the result of multiplying \code{y} on the left
      by \eqn{P_{1}' Q}{P1' * Q}.}
    \item{\code{solve}}{\code{signature(a = "sparseQR", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\details{
  For a sparse \eqn{m \times n}{m x n} (\dQuote{long}: \eqn{m \ge n}{m >= n})
  rectangular matrix \eqn{A}, the sparse QR decomposition is either 
  \cr of the form \eqn{P A = Q R} with a (row)
  permutation matrix \eqn{P}, (encoded in the \code{p} slot of the
  result) if the \code{q} slot is of length 0,
  \cr or of the form \eqn{P A P* = Q R} with an extra (column) permutation
  matrix \eqn{P*} (encoded in the \code{q} slot).
  Note that the row permutation \eqn{P A} in \R is simply \code{A[p+1, ]}
  where \code{p} is the \code{p}-slot, a 0-based permutation of
  \code{1:m} applied to the rows of the original matrix.

  If the \code{q} slot has length \code{n} it is a 0-based permutation
  of \code{1:n} applied to the columns of the original matrix to reduce
  the amount of \dQuote{fill-in} in the matrix \eqn{R}, and
  \eqn{A P*} in \R is simply \code{A[ , q+1]}.

  \eqn{R} is an \eqn{m\times n}{m by n} matrix that is zero below the
  main diagonal, i.e., upper triangular (\eqn{m\times m}{m by n}) with
  \eqn{m-n} extra zero rows.

  The matrix \eqn{Q} is a \dQuote{virtual matrix}.  It is the product of
  \eqn{n} Householder transformations.  The information to generate
  these Householder transformations is stored in the \code{V} and
  \code{beta} slots.
  \cr
  Note however that \code{qr.Q()} returns the row permuted matrix
  \eqn{Q* := P^{-1}Q = P'Q}{Q* := P^(-1) Q = P'Q} as permutation matrices are
  orthogonal; and \eqn{Q*} is orthogonal itself because \eqn{Q} and \eqn{P} are.
  This is useful because then, as in the dense matrix and \pkg{base} \R
  matrix \code{\link{qr}} case, we have the mathematical identity
  \deqn{P A = Q* R,} in \R as \preformatted{            A[p+1,] == qr.Q(*) \%*\% R .}

  The \code{sparseQR} methods for the \code{qr.*} functions return
  objects of class \code{\linkS4class{dgeMatrix}}.
  Results from \code{qr.coef},
  \code{qr.resid} and \code{qr.fitted} (when \code{k == ncol(R)}) are
  well-defined and should match those from the corresponding dense matrix
  calculations.  However, because the matrix \code{Q} is not uniquely
  defined, the results of \code{qr.qy} and \code{qr.qty} do not
  necessarily match those from the corresponding dense matrix
  calculations.

  Also, the results of \code{qr.qy} and \code{qr.qty} apply to the
  permuted column order when the \code{q} slot has length \code{n}.
}
\seealso{
  Class \code{\linkS4class{dgCMatrix}}.
  
  Generic function \code{\link[base]{qr}} from \pkg{base},
  whose default method \code{qr.default} \dQuote{defines}
  the S3 class \code{qr} of dense QR factorizations.

  \code{\link{qr-methods}} for methods defined in \pkg{Matrix}.
  
  The many auxiliary functions for QR factorizations:
  \code{\link{qr.Q}}, \code{\link{qr.R}}, \code{\link{qr.X}},
  \code{\link{qr.coef}}, \code{\link{qr.fitted}},
  \code{\link{qr.qty}}, \code{\link{qr.qy}},
  \code{\link{qr.resid}}, and \code{\link{qr.solve}}.
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
