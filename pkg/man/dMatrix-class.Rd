\name{dMatrix-class}
\docType{class}
\alias{dMatrix-class}
\alias{show,dMatrix-method}
\alias{coerce,matrix,dMatrix-method}
% Group methods
\alias{Arith,dMatrix,dMatrix-method}
\alias{Arith,dMatrix,numeric-method}
\alias{Arith,numeric,dMatrix-method}
\alias{Math,dMatrix-method}
\alias{Math2,dMatrix,numeric-method}% Math2 = round + signif, but
% for silly reasons, need these 2 as well:
\alias{round,dMatrix,numeric-method}
\alias{signif,dMatrix,numeric-method}
%
\title{(Virtual) Class "dMatrix" of "double" Matrices}
\description{
  The \code{dMatrix} class is a virtual class contained by all actual
  classes in the \pkg{dMatrix} package.
}
%% FIXME: It's no longer virtual in the strict sense since it *has* slots
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
  The are (very simplistic) group methods (see, e.g., \code{\link{Arith}})
  \describe{
    \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "dMatrix")}: ... }
    \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "dMatrix")}: ... }
    \item{Math}{\code{signature(x = "dMatrix")}: ... }
    \item{Math2}{\code{signature(x = "dMatrix", digits = "numeric")}: ...}
  }
}
%\references{}
% Martin + Doug\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{dgCMatrix-class}}, and
  \code{\link{Matrix-class}}.
}
\examples{
 M <- new("dMatrix")
 (M <- new("dMatrix", Dim = 1:2, Dimnames = list("r1", c("c1","c2"))))
 str(M)
}
\keyword{classes}
\keyword{algebra}
