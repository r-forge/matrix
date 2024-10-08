\name{model.Matrix}
\alias{model.Matrix}
%% Cut and paste a lot from the "standard"  model.matrix() help page,
%% <--> ~/R/D/r-devel/R/src/library/stats/man/model.matrix.Rd
\title{Construct Possibly Sparse Design or Model Matrices}
\description{
  \code{model.Matrix} creates design matrix, very much like the
  standard \R function \code{\link{model.matrix}}, however returning a
  dense or sparse object of class \code{\linkS4class{modelMatrix}}.
}
\usage{
model.Matrix(object, data = environment(object),
             contrasts.arg = NULL, xlev = NULL,
             sparse = FALSE, drop.unused.levels = FALSE, \dots)
}
\arguments{
  \item{object}{an object of an appropriate class.  For the default
    method, a model \link{formula} or a \code{\link{terms}} object.}
  \item{data}{a data frame created with \code{\link{model.frame}}.  If
    another sort of object, \code{model.frame} is called first.}
  \item{contrasts.arg}{A list, whose entries are values (numeric
    matrices or character strings naming functions) to be used
    as replacement values for the \code{\link{contrasts}}
    replacement function and whose names are the names of
    columns of \code{data} containing \code{\link{factor}}s.}
  \item{xlev}{to be used as argument of \code{\link{model.frame}} if
    \code{data} has no \code{"terms"} attribute.}
  \item{sparse}{logical indicating if the result should be sparse
    (of class \code{\linkS4class{sparseModelMatrix}}), using
    \code{\link[Matrix]{sparse.model.matrix}()} (package \pkg{Matrix}).}
  \item{drop.unused.levels}{used only when \code{sparse} is TRUE: Should
    factors have unused levels dropped?
    (This used to be true, \emph{implicitly} in the first versions up to
    July 2010; the default has been changed for compatibility with
    \R's standard (dense) \code{\link{model.matrix}()}.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{model.Matrix()} is a simple wrapper either (\code{sparse = FALSE})
  around the traditional \code{\link{model.matrix}()} returning a
  \code{"\linkS4class{ddenseModelMatrix}"}, or (\code{sparse = TRUE})
  around \code{\link[Matrix]{sparse.model.matrix}()}, returning a
  \code{"\linkS4class{dsparseModelMatrix}"} object.

  \code{model.Matrix} creates a design matrix from the description given
  in \code{terms(object)}, using the data in \code{data} which must
  supply variables with the same names as would be created by a call to
  \code{model.frame(object)} or, more precisely, by evaluating
  \code{attr(terms(object), "variables")}.

  For more details, see \code{\link{model.matrix}}.
}
\value{
 an object inheriting from class \code{\linkS4class{modelMatrix}}, by
 default, \code{\linkS4class{ddenseModelMatrix}}.
}
% \author{Martin Maechler}
\seealso{
  \code{\link[stats]{model.matrix}}, and
  \code{\link[Matrix]{sparse.model.matrix}} from package \CRANpkg{Matrix}.
}
\examples{
data(CO2, package="datasets")
class(sm <- model.Matrix(~ 0+Type*Treatment, data=CO2, sparse=TRUE))
class(dm <- model.Matrix(~ 0+Type*Treatment, data=CO2, sparse=FALSE))
stopifnot(dim(sm) == c(84,4), dim(sm) == dim(dm), all(sm == dm))
}
\keyword{models}
