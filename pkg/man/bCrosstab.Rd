\name{bCrosstab}
\alias{bCrosstab}
\title{Create pairwise crosstabulation}
\description{
  Create the pairwise crosstabulation of two or more factors and store
  it as an object of class \code{"bCrosstab"} which consists of a
  list of sparse, blocked matrices.
}
\usage{
bCrosstab(flist, reorder = TRUE)
}
\arguments{
  \item{flist}{A list of two or more factors.  All factors must have the
    same length.}
  \item{reorder}{Logical flag determining if reordering of the factors
    is allowed.}
}
%\details{}
\value{
  A list.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{bbCrosstab-class}}}
%\examples{}
\keyword{manip}
\keyword{algebra}
