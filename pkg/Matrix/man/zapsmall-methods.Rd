\name{zapsmall-methods}
%%
%% \alias{zapsmall}
\alias{zapsmall-methods}
\alias{zapsmall,Matrix-method}
\alias{zapsmall,sparseVector-method}
%%
\title{Rounding of Numbers Zapping Small Ones to Zero}
\description{
  S4 methods defined in package \pkg{Matrix} for S4 generic function
  \code{zapsmall}, generalizing \code{base::\link{zapsmall}}.

  The S4 generic function definition is obtained from package
  \pkg{methods} through the \code{\link{implicitGeneric}} mechanism.
}
\usage{
\S4method{zapsmall}{Matrix}(x, digits = getOption("digits"),
         mFUN = function(x, ina) max(abs(x[!ina])),
         min.d = 0L, \dots)
\S4method{zapsmall}{sparseVector}(x, digits = getOption("digits"),
         mFUN = function(x, ina) max(abs(x[!ina])),
         min.d = 0L, \dots)
}
\arguments{
  \item{x}{a \code{\linkS4class{Matrix}} object or
    \code{\linkS4class{sparseVector}} object.  Pattern, logical, and
    integer data are converted to double precision floating-point.}
  \item{digits, mFUN, min.d}{additional arguments passed to
    \code{base::\link{zapsmall}}.}
  \item{\dots}{optional arguments not yet used.}
}
\details{
  Methods preserve the nonzero pattern of sparse format \code{x}: zeros
  introduced by rounding are nonstructural.
}
\value{
  An object representing the result of rounding \code{x}.  The class is
  a superclass of the class of \code{x} unless \code{x} has pattern,
  logical, or integer data.
}
\seealso{
  \code{\link{drop0}} for deleting small (on an absolute scale) entries.
}
\examples{
x <- new("dgCMatrix", Dim = c(3L, 2L),
         p = c(0L, 3L, 4L), i = c(0L, 1L, 2L, 1L),
         x = c(0, 1, -51, 1e+6)) # maximum absolute value is 10^6
x
zapsmall(x, digits = 4L, min.d = -Inf) # rounds to nearest 10^(6-4)
zapsmall(x, digits = 3L, min.d = -Inf) # rounds to nearest 10^(6-3)
}
\docType{methods}
\keyword{array}
\keyword{manip}
\keyword{methods}
