\name{unpack}
\title{Full Storage Representation of Packed Matrices}
\usage{
unpack(x, \dots)
}
\alias{unpack}
\description{
  Expands decompositions stored in compact form into matrix factors.
}
\arguments{
  \item{x}{
    a matrix stored in packed form.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  A \code{Matrix} object containing the full-storage representation of \code{x}.
}
\details{
  This is a generic function with special methods for different types
  of packed matrices.
  Use \code{\link{methods}("unpack")} to list all the methods for the
  \code{unpack} generic.
}
\examples{
showMethods("unpack")
cp4 <- chol(Hilbert(4))
tp4 <- as(cp4,"dtpMatrix")
str(tp4)
(unpack(tp4))

\dontrun{
## once we have Diagonal() :
unpack(Diagonal( 1:3))
}
}
\keyword{array}
\keyword{algebra}
