\name{sparseMatrix-class}
\title{Virtual Class \sQuote{sparseMatrix} of Sparse Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{sparseMatrix-class}
%
\description{
  \code{sparseMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing sparse format matrices.
  Formally, a sparse format is a storage format in which the space
  required to represent an \eqn{m \times n}{m-by-n} matrix with \eqn{p}
  nonzero entries is \eqn{O(\max\{m, n, p\})}{O(max(m, n, p))}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
  }
}
\section{Subclasses}{
  Package \pkg{Matrix} defines two virtual subclasses of
  \code{denseMatrix}: \code{\linkS4class{unpackedMatrix}} and
  \code{\linkS4class{packedMatrix}}.  These correspond to the
  conventional rectangular and packed storage formats used by LAPACK:
  \url{https://netlib.org/lapack/lug/node121.html}.  See those help
  topics for implementation details and (recursively) further
  subclasses.
}
\seealso{
  \dQuote{Complementary} class \code{\linkS4class{denseMatrix}}
  representing dense (meaning \emph{not} sparse) format matrices.
  Function \code{\link{sparseMatrix}}, a widely used (but not the only)
  constructor for \code{sparseMatrix} objects.  Other options for
  construction, not including explicit coercion:
  \code{\link{Matrix}(sparse = TRUE)},
  \code{\link{xtabs}(sparse = TRUE)},
  \code{\link{sparse.model.matrix}()}, \ldots
}
\examples{
showClass("sparseMatrix") ## and look at the help() of its subclasses
M <- Matrix(0, 10000, 100)
M[1,1] <- M[2,3] <- 3.14
M  ## show(.) method suppresses printing of the majority of rows

data(CAex, package = "Matrix")
dim(CAex) # 72 x 72 matrix
determinant(CAex) # works via sparse lu(.)

## factor -> t( <sparse design matrix> ) :
(fact <- gl(5, 3, 30, labels = LETTERS[1:5]))
(Xt <- as(fact, "sparseMatrix"))  # indicator rows

## missing values --> all-0 columns:
f.mis <- fact
i.mis <- c(3:5, 17)
is.na(f.mis) <- i.mis
Xt != (X. <- as(f.mis, "sparseMatrix")) # differ only in columns 3:5,17
stopifnot(all(X.[,i.mis] == 0), all(Xt[,-i.mis] == X.[,-i.mis]))
}
