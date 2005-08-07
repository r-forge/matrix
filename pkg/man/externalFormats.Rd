\name{externalFormats}
\alias{readHB}
\alias{readMM}
\title{Read and write external matrix formats}
\description{
    Read matrices stored in the Harwell-Boeing or MatrixMarket formats.
}
\usage{
readHB(filename)
readMM(filename)
}
\arguments{
  \item{filename}{the name of the file to read, as a character scalar.
    The names of files storing matrices in the Harwell-Boeing format
    usually end in \code{".rua"} or \code{".rsa"}.  Those storing
    matrices in the MatrixMarket format usually end in {".mtx"}.}
}
\value{
  An object that inherits from the \code{"Matrix"} class.
}
\references{
  \url{http://math.nist.gov/MatrixMarket}
}
\keyword{IO}
\keyword{algebra}

