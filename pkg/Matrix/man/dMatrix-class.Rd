\name{dMatrix-class}
\title{Virtual Class "dMatrix" of Double Precision Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{dMatrix-class}
\alias{ddenseMatrix-class}
\alias{dsparseMatrix-class}
%
\alias{coerce,matrix,dMatrix-method}
\alias{coerce,vector,dMatrix-method}
%
\alias{coerce,matrix,ddenseMatrix-method}
\alias{coerce,vector,ddenseMatrix-method}
%
\alias{coerce,matrix,dsparseMatrix-method}
\alias{coerce,vector,dsparseMatrix-method}
%
\description{
  The \code{dMatrix} class is a virtual class contained by all actual
  classes of numeric matrices in the \pkg{Matrix} package.  Similarly,
  all the actual classes of logical matrices inherit from the
  \code{lMatrix} class.
}
%\section{Objects from the Class}{A virtual Class: No objects may be
%  created from it.
%}
\section{Slots}{
  Common to \emph{all} matrix object in the package:
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Methods}{
  There are (relatively simple) group methods (see, e.g., \code{\link{Arith}})
  \describe{
    \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "dMatrix")}: ... }
    \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "dMatrix")}: ... }
    \item{Math}{\code{signature(x = "dMatrix")}: ... }
    \item{Math2}{\code{signature(x = "dMatrix", digits = "numeric")}:
      this group contains \code{\link{round}()} and \code{\link{signif}()}.}
    \item{Compare}{\code{signature(e1 = "numeric", e2 = "dMatrix")}: ... }
    \item{Compare}{\code{signature(e1 = "dMatrix", e2 = "numeric")}: ... }
    \item{Compare}{\code{signature(e1 = "dMatrix", e2 = "dMatrix")}: ... }
    \item{Summary}{\code{signature(x = "dMatrix")}: The \code{"Summary"}
      group contains the seven functions
      \code{\link{max}()}, \code{\link{min}()}, \code{\link{range}()},
      \code{\link{prod}()}, \code{\link{sum}()},
      \code{\link{any}()}, and \code{\link{all}()}.}
  }
  The following methods are also defined for all double matrices:
  \describe{
    \item{expm}{\code{signature(x = "dMatrix")}: computes the
      \emph{\dQuote{Matrix Exponential}}, see \code{\link{expm}}.}
  }
  The following methods are defined for all logical matrices:
  \describe{
    \item{which}{\code{signature(x = "lsparseMatrix")} and many other
      subclasses of \code{"lMatrix"}: as the \pkg{base} function
      \code{\link{which}(x, arr.ind)} returns the indices of the
      \code{\link{TRUE}} entries in \code{x}; if \code{arr.ind} is true,
      as a 2-column matrix of row and column indices. Since \pkg{Matrix}
      version 1.2-9, if \code{useNames} is true, as by default, with
      \code{\link{dimnames}}, the same as \code{base::which}.}
  }
}
%\references{}
% Martin + Doug\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  The nonzero-pattern matrix class \code{\linkS4class{nMatrix}}, which
  can be used to store non-\code{\link{NA}} \code{\link{logical}}
  matrices even more compactly.

  The numeric matrix classes \code{\linkS4class{dgeMatrix}},
  \code{\linkS4class{dgCMatrix}}, and \code{\linkS4class{Matrix}}.

  \code{\link{drop0}(x, tol=1e-10)} is sometimes preferable to (and
  more efficient than) \code{zapsmall(x, digits=10)}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
 showClass("dMatrix")

 set.seed(101)
 round(Matrix(rnorm(28), 4,7), 2)
 M <- Matrix(rlnorm(56, sd=10), 4,14)
 (M. <- zapsmall(M))
 table(as.logical(M. == 0))

showClass("ddenseMatrix")

showMethods(class = "ddenseMatrix", where = "package:Matrix")

showClass("dsparseMatrix")
}
