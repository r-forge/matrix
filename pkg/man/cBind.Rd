\name{cBind}
\alias{cBind}
\alias{rBind}
\alias{cbind}
\alias{rbind}
%
\title{Versions of 'cbind' and 'rbind' recursively built on cbind2/rbind2}
\description{
  .................
  ..................

  FIXME

  ..................
}
\usage{
cbind(..., deparse.level = 1)
rbind(..., deparse.level = 1)

cBind(..., deparse.level = 1)
rBind(..., deparse.level = 1)
}
\arguments{
  \item{\dots}{matrix-like \R objects to be bound together, see
    \code{\link{cbind}} and \code{\link{rbind}}.}
  \item{deparse.level}{integer determining under which circumstances
    column and row names are built from the actual arguments'
    \sQuote{expression}, see \code{\link{cbind}}.}
}
\details{
  The implementation of these is \emph{recursive}, calling
  \code{\link[methods]{cbind2}} or
  \code{\link[methods]{rbind2}} respectively, where these have methods
  defined and so should dispatch appropriately.
}
\value{
  typically a \sQuote{matrix-like} object of a similar
  \code{\link{class}} as the first argument in \code{\dots}.
}
% \references{ ~put references to the literature/web site here ~ }
% \author{ ~~who you are~~ }
% \note{ ~~further notes~~

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{cbind2}}, \code{\link{cbind}},
  \code{\link[methods]{Methods}}.
}
\examples{
%%%%%%% TODO
}
\keyword{array}
\keyword{manip}

