\name{bbCrosstab-class}
\docType{class}
\alias{bbCrosstab-class}
\alias{image,bbCrosstab-method}
\title{Pairwise crosstabulations}
\description{An \code{bbCrosstab} object is a pairwise crosstabulation
  of two or more factors expressed as list of block/block compressed sparse matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bbCrosstab", ...)}
or by a call to the constructor function \code{bbCrosstab}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{named list of \code{cscBlocked} objects.}
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
  }
}
\section{Extends}{
Class \code{"list"}, directly.
}
\section{Methods}{
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\code{\link{bbCrosstab}}}
}
%\examples{}
\keyword{classes}
