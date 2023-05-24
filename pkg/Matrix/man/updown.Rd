\name{updown-methods}
\title{Up- and Down-Dating a Cholesky Decomposition}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{updown}
\alias{updown-methods}
%
\alias{updown,character,Matrix,CHMfactor-method}
\alias{updown,character,matrix,CHMfactor-method}
\alias{updown,logical,Matrix,CHMfactor-method}
\alias{updown,logical,dgCMatrix,CHMfactor-method}
\alias{updown,logical,dsCMatrix,CHMfactor-method}
\alias{updown,logical,dtCMatrix,CHMfactor-method}
\alias{updown,logical,matrix,CHMfactor-method}
%
\description{
  Compute the up- or down-dated Cholesky decomposition
}
\usage{
updown(update, C, L)
}
\arguments{
  \item{update}{logical (\code{TRUE} or \code{FALSE}) or \code{"+"} or
    \code{"-"} indicating if an up- or a down-date is to be computed.}
  \item{C}{any \R object, coercable to a sparse matrix (i.e., of
    subclass of \code{\linkS4class{sparseMatrix}}).}
  \item{L}{a Cholesky factor, specifically, of class
    \code{"\linkS4class{CHMfactor}"}.}
}

% \details{
% }

\value{
  an updated Cholesky factor, of the same dimension as \code{L}.
  Typically of class \code{"\linkS4class{dCHMsimpl}"} (a sub class of
  \code{"\linkS4class{CHMfactor}"}).
}
\references{
  CHOLMOD manual, currently beginning of chapter~18. ... %% FIXME
}
\author{Contributed by Nicholas Nagle, University of Tennessee, Knoxville, USA}
\section{Methods}{
  \describe{
    \item{\code{signature(update = "character", C = "mMatrix", L = "CHMfactor")}}{..}
    \item{\code{signature(update = "logical", C = "mMatrix", L = "CHMfactor")}}{ .. }
  }
}
\seealso{
 \code{\link{Cholesky}},
}
\examples{
m <- sparseMatrix(i = c(3, 1, 3:2, 2:1), p = c(0:2, 4, 4, 6), x = 1:6,
                  dimnames = list(LETTERS[1:3], letters[1:5]))
uc0 <- Cholesky(A <- crossprod(m) + Diagonal(5))
uc1 <- updown("+", Diagonal(5, 1), uc0)
uc2 <- updown("-", Diagonal(5, 1), uc1)
stopifnot(all.equal(uc0, uc2))

\dontshow{
if(FALSE) {
## Hmm: this loses positive definiteness:
uc2 <- updown("-", Diagonal(5, 2), uc0)
image(show(as(uc0, "CsparseMatrix")))
image(show(as(uc2, "CsparseMatrix"))) # severely negative entries
}
}
}
