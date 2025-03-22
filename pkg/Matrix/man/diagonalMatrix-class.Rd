\name{diagonalMatrix-class}
\title{Virtual Class \sQuote{diagonalMatrix} of Diagonal Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{diagonalMatrix-class}
%
\alias{ndiMatrix-class}
\alias{ldiMatrix-class}
\alias{idiMatrix-class}
\alias{ddiMatrix-class}
\alias{zdiMatrix-class}
%
\description{
  \code{diagonalMatrix} is a virtual subclass of
  \code{\linkS4class{triangularMatrix}} and of
  \code{\linkS4class{indexMatrix}} representing diagonal matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.  \code{Dim[1]} and \code{Dim[2]} must
      be equal.}
  }
  All direct, nonvirtual subclasses of \code{diagonalMatrix} defined in
  package \pkg{Matrix} have two additional slots:
  \describe{
    \item{\code{diag}}{a character string, either \code{"N"} or
      \code{"U"}.  \code{"N"} indicates that \code{x} stores the
      diagonal entries of the represented matrix.  \code{"U"} indicates
      that the represented matrix is unit diagonal and that \code{x} is
      empty.}
    \item{\code{x}}{a vector of length \code{n=Dim[1]=Dim[2]} or 0,
      conditionally (depending on \code{diag}) storing the diagonal
      entries of the represented matrix.}
  }
}
\section{Subclasses}{
  The direct, nonvirtual subclasses of \code{diagonalMatrix} defined in
  package \pkg{Matrix} are listed below with links to virtual
  superclasses defining their data type.
  \tabular{rr}{
    Nonvirtual subclass \tab Data type\cr
    \code{ndiMatrix} \tab \code{\linkS4class{nMatrix}}\cr
    \code{ldiMatrix} \tab \code{\linkS4class{lMatrix}}\cr
    \code{idiMatrix} \tab \code{\linkS4class{iMatrix}}\cr
    \code{ddiMatrix} \tab \code{\linkS4class{dMatrix}}\cr
    \code{zdiMatrix} \tab \code{\linkS4class{zMatrix}}
  }
  Class \code{ndiMatrix} was not defined until \pkg{Matrix} version
  1.6-2.
}
\seealso{
  Related classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{posdefMatrix}},
  \code{\linkS4class{triangularMatrix}}, and
  \code{\linkS4class{indexMatrix}} representing
  general (unstructured, possibly nonsquare), Hermitian or symmetric,
  positive semidefinite, triangular, and index matrices.  Generic
  function \code{\linkS4methods{isDiagonal}} for testing for testing
  if a matrix is diagonal.  Generic function
  \code{\linkS4methods{forceDiagonal}} for inducing diagonal structure.
  Function \code{\link{Diagonal}}, a constructor for
  \code{diagonalMatrix} objects.
}
\examples{
I5 <- Diagonal(5)
D5 <- Diagonal(x = 10*(1:5))
## trivial (but explicitly defined) methods:
stopifnot(identical(crossprod(I5), I5),
          identical(tcrossprod(I5), I5),
          identical(crossprod(I5, D5), D5),
          identical(tcrossprod(D5, I5), D5),
          identical(solve(D5), solve(D5, I5)),
          all.equal(D5, solve(solve(D5)), tolerance = 1e-12)
          )
solve(D5)# efficient as is diagonal

# an unusual way to construct a band matrix:
rbind2(cbind2(I5, D5),
       cbind2(D5, I5))

%% moved from ./ldiMatrix-class.Rd :
(lM <- Diagonal(x = c(TRUE,FALSE,FALSE)))
str(lM)#> gory details (slots)

crossprod(lM) # numeric
(nM <- as(lM, "nMatrix"))
crossprod(nM) # pattern sparse

%% moved from ./ddiMatrix-class.Rd :
(d2 <- Diagonal(x = c(10,1)))
str(d2)
## slightly larger in internal size:
str(as(d2, "sparseMatrix"))

M <- Matrix(cbind(1,2:4))
M \%*\% d2 #> `fast' multiplication

chol(d2) # trivial
stopifnot(is(cd2 <- chol(d2), "ddiMatrix"),
          all.equal(cd2@x, c(sqrt(10),1)))
}
