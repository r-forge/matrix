\name{Matrix-package}
\title{\R{} Package \pkg{Matrix}}
\docType{package}
\keyword{package}
\alias{Matrix-package}
\alias{Matrix.Version}
\description{TODO.}
\usage{
Matrix.Version()
}
%% \details{}
%% \author{}
\examples{
(Mv <- Matrix.Version())
stopifnot(identical(Mv[["package"]], utils::packageVersion("Matrix")))
}
