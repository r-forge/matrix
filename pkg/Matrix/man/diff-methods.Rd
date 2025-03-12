\name{diff-methods}
\title{Difference a Vector}
%
\docType{methods}
\keyword{arith}
\keyword{methods}
%
\alias{diff}
\alias{diff-methods}
%
\alias{diff,denseMatrix-method}
\alias{diff,sparseMatrix-method}
\alias{diff,sparseVector-method}
%
\description{TODO.}
\usage{
diff(x, \dots)
\S4method{diff}{sparseVector}(x, lag = 1L, differences = 1L, \dots)
}
\arguments{
\item{x}{.}
\item{lag}{.}
\item{differences}{.}
\item{\dots}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

%% \item{diff}{\code{signature(x = "Matrix")}: As \code{\link{diff}()}
%%   for traditional matrices, i.e., applying \code{diff()} to each
%%   column.}

\examples{
(M <- Matrix(c(0, 1, 0, 0), 6, 4))
diff(M)
}
