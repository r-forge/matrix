\name{Cholesky-class}
\title{Cholesky and Bunch-Kaufman Decompositions}
\docType{class}
%
\alias{Cholesky-class}
\alias{pCholesky-class}
\alias{BunchKaufman-class}
\alias{pBunchKaufman-class}
%
\alias{coerce,BunchKaufman,Matrix-method}
\alias{coerce,pBunchKaufman,Matrix-method}
\alias{show,BunchKaufman-method}
\alias{show,pBunchKaufman-method}
\alias{show,Cholesky-method}
\alias{show,pCholesky-method}
%
\description{The \code{"Cholesky"} class is the class of Cholesky
  decompositions of positive-semidefinite, real dense matrices.  The
  \code{"BunchKaufman"} class is the class of Bunch-Kaufman
  decompositions of symmetric, real matrices.  The \code{"pCholesky"}
  and \code{"pBunchKaufman"} classes are their \emph{\bold{p}acked}
  storage versions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Cholesky",
    ...)} or \code{new("BunchKaufman", ...)}, etc,
  or rather by calls of the form \code{\link{chol}(pm)} or
  \code{\link{BunchKaufman}(pm)} where \code{pm} inherits from
  the \code{"\linkS4class{dpoMatrix}"} or \code{"\linkS4class{dsyMatrix}"}
  class or as a side-effect of other functions
  applied to \code{"dpoMatrix"} objects (see \code{\linkS4class{dpoMatrix}}).
}
\note{
  \enumerate{
    \item Objects of class \code{"Cholesky"} typically stem from
      \code{\link{chol}(D)}, applied to a \emph{dense} matrix \code{D}.

      On the other hand, the \emph{function} \code{\link{Cholesky}(S)}
      applies to a \emph{sparse} matrix \code{S}, and results in objects
      inheriting from class \code{\linkS4class{CHMfactor}}.

    \item For traditional matrices \code{m}, \code{chol(m)} is a
      traditional matrix as well, triangular, but simply an \eqn{n\times
	n}{n * n} numeric \code{\link{matrix}}.  Hence, for compatibility,
      the \code{"Cholesky"} class extends a triangular Matrix class,
      namely "dtrMatrix", and, analogously, \code{"pCholesky"} is
      defined as an extension of \code{"dtpMatrix"}.

      Consequently, \code{\link{determinant}(R)} for \code{R <- chol(A)}
      returns the determinant of \code{R}, not of \code{A}.  This is in
      contrast to class \code{\linkS4class{CHMfactor}} objects \code{C},
      where \code{determinant(C)} gives the determinant of the
      \emph{original} matrix \code{A}, for \code{C <- Cholesky(A)},
      see also the \code{determinant} method documentation on the class
      \code{\linkS4class{CHMfactor}} page.

    \item Although the internal representation of class
      \code{"BunchKaufman"} is triangular, similar to \code{"Cholesky"},
      the class does not extend \code{"dtrMatrix"}, because, unlike
      \code{"Cholesky"}, what is represented is multiple matrix factors
      in a condensed format.  To obtain a \code{"dtrMatrix"}, one can
      use \code{as(., "Matrix")}.
  }
}
%% FIXME: all _incredibly_ repetitive ...
\section{Slots}{
  A Cholesky decomposition extends class
  \code{\linkS4class{MatrixFactorization}} but is basically a triangular
  matrix extending the \code{"\linkS4class{dtrMatrix}"} class.
  \describe{
    \item{\code{uplo}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{diag}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{x}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{Dim}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{Dimnames}:}{inherited from the \code{"dtrMatrix"} class.}
  }
  The packed version extends the \code{"dtpMatrix"} class.
  
  A Bunch-Kaufman decomposition has all of the slots of a
  \code{"dtrMatrix"} and \code{"dtpMatrix"} without extending
  either class (or \code{"Matrix"} whatsoever).  It also has
  a \code{perm} slot representing a permutation matrix.
}
\section{Extends}{
%% FIXME: only true for Cholesky ...
  Class \code{"MatrixFactorization"} and \code{"dtrMatrix"}, directly.
  Class \code{"Matrix"}, by class \code{"dtrMatrix"}.
}
\section{Methods}{
  As \code{"Cholesky"} extends \code{"dtrMatrix"}, it is treated
  \emph{directly} as a (triangular) matrix in computations.  That is
  not the case for \code{"BunchKaufman"}, for which (currently)
  there is only \code{show()} and the coercion \code{as(., "Matrix")}.
}
\seealso{
  Classes \code{\linkS4class{dtrMatrix}}, \code{\linkS4class{dpoMatrix}};
  function \code{\link{chol}}.

  Function \code{\link{Cholesky}} resulting in class
  \code{\linkS4class{CHMfactor}} objects, \emph{not} class "Cholesky"
  ones, see the section \sQuote{Note}.
}
\examples{
(sm <- pack(Matrix(diag(5) + 1))) # dspMatrix
signif(csm <- chol(sm), 4)

(pm <- crossprod(Matrix(rnorm(18), nrow = 6, ncol = 3)))
(ch <- chol(pm))
if (toupper(ch@uplo) == "U") # which is TRUE
   crossprod(ch)
stopifnot(all.equal(as(crossprod(ch), "matrix"),
                    as(pm, "matrix"), tolerance=1e-14))
}
\keyword{classes}
\keyword{algebra}
