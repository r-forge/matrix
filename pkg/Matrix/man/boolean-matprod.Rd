\name{\%&\%-methods}
\title{Boolean Arithmetic Matrix Products: \code{\%&\%} and Methods}
\docType{methods}
\alias{\%&\%}
\alias{\%&\%-methods}
\alias{\%&\%,ANY,ANY-method}
\alias{\%&\%,ANY,matrix-method}
\alias{\%&\%,ANY,Matrix-method}
\alias{\%&\%,matrix,ANY-method}
\alias{\%&\%,Matrix,ANY-method}
\alias{\%&\%,matrix,matrix-method}
\alias{\%&\%,Matrix,Matrix-method}
\alias{\%&\%,mMatrix,mMatrix-method}
\alias{\%&\%,mMatrix,nMatrix-method}
\alias{\%&\%,nMatrix,mMatrix-method}
\alias{\%&\%,nMatrix,nMatrix-method}
\alias{\%&\%,nMatrix,nsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nMatrix-method}
\alias{\%&\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,nCsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,nCsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,geMatrix-method}
\alias{\%&\%,geMatrix,diagonalMatrix-method}
\alias{\%&\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,sparseMatrix-method}
\alias{\%&\%,sparseMatrix,diagonalMatrix-method}

\description{
  For boolean or \dQuote{patter\bold{n}} matrices, i.e., \R object of
  class \code{\linkS4class{nMatrix}}, it is natural to allow matrix
  products using boolean instead of numerical arithmetic.

  In package \pkg{Matrix}, we use the (\dQuote{infix}) function aka
  binary operator \code{\%&\%} for this and provide methods for all our
  matrices and the traditional \R matrices (see \code{\link{matrix}}).
}
\section{Methods}{
  We provide methods for both the \dQuote{traditional} (\R base) matrices
  and conceptually all matrices in package \pkg{Matrix}.
  \describe{
    \item{\code{signature(x = "ANY", y = "ANY")}}{
    }

    \item{\code{signature(x = "ANY", y = "matrix")}}{
    }

    \item{\code{signature(x = "ANY", y = "Matrix")}}{
    }

    \item{\code{signature(x = "matrix", y = "ANY")}}{
    }

    \item{\code{signature(x = "Matrix", y = "ANY")}}{
    }

    \item{\code{signature(x = "matrix", y = "matrix")}}{
    }

    \item{\code{signature(x = "Matrix", y = "Matrix")}}{
    }

    \item{\code{signature(x = "Matrix", y = "nMatrix")}}{
    }

    \item{\code{signature(x = "nMatrix", y = "Matrix")}}{
    }

    \item{\code{signature(x = "nMatrix", y = "nMatrix")}}{
    }

    \item{\code{signature(x = "nMatrix", y = "nsparseMatrix")}}{
    }

    \item{\code{signature(x = "nsparseMatrix", y = "nMatrix")}}{
    }

    \item{\code{signature(x = "nsparseMatrix", y = "nsparseMatrix")}}{
    }
  }% {describe}
}% {Methods}
\value{
  a pattern matrix, i.e., inheriting from \code{"\linkS4class{nMatrix}"},
  or an \code{"\linkS4class{ldiMatrix}"} in case of a diagonal matrix.
}
\examples{
set.seed(7)
L <- Matrix(rnorm(20) > 1,    4,5)
(N <- as(L, "nMatrix"))
D <- Matrix(round(rnorm(30)), 5,6) # -> values in -1:1 (for this seed)
L \%&\% D
stopifnot(identical(L \%&\% D, N \%&\% D),
          all(L \%&\% D == as((L \%*\% abs(D)) > 0, "sparseMatrix")))
crossprod(N)     # -> sparse patter'n' (TRUE/FALSE : boolean arithmetic)
crossprod(N  +0) # -> numeric Matrix (with same "pattern")
stopifnot(all(crossprod(N) == t(N) \%&\% N))
crossprod(D, boolArith =  TRUE) # <- currently with warning "... not yet implemented":
crossprod(L, boolArith =  TRUE) # pattern: "nsCMatrix"
crossprod(L, boolArith = FALSE) # numeric: "dsCMatrix"
}
\keyword{methods}
