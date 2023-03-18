\name{FastCoercions}
\title{}
%
\alias{FastCoercions}
%
\alias{.CR2RC}
\alias{.CR2T}
\alias{.M2diag}
\alias{.M2sym}
\alias{.M2tri}
\alias{.T2CR}
\alias{.dense2g}
\alias{.dense2kind}
\alias{.dense2m}
\alias{.dense2sparse}
\alias{.diag2kind}
\alias{.diag2m}
\alias{.diag2sparse}
\alias{.m2dense}
\alias{.m2sparse}
\alias{.sparse2dense}
\alias{.sparse2g}
\alias{.sparse2kind}
\alias{.sparse2m}
\alias{.tCR2RC}
%
\description{
}
\usage{
.M2tri(from, ...)
.M2sym(from, ...)
.M2diag(from, check = TRUE)

.m2dense(from, code, uplo = "U", diag = "N")
.m2sparse(from, code, uplo = "U", diag = "N")

.dense2m(from, ndense = FALSE)
.sparse2m(from)
.diag2m(from)

.dense2kind(from, kind)
.sparse2kind(from, kind, drop0 = FALSE)
.diag2kind(from, kind)

.dense2g(from, kind = ".")
.sparse2g(from)

.sparse2dense(from, packed = FALSE)

.dense2sparse(from, code)

.diag2sparse(from, code, uplo = "U", drop0 = TRUE)

.CR2T(from)
.T2CR(from, Csparse = TRUE)

.CR2RC(from)
.tCR2RC(from)
}
\arguments{
  \item{from}{}
  \item{check}{}
  \item{code}{}
  \item{kind}{}
  \item{uplo}{}
  \item{diag}{}
  \item{drop0}{}
  \item{packed}{}
  \item{ndense}{}
  \item{Csparse}{}
  \item{\dots}{}
}
\value{
}
\details{
}
\seealso{
}
\keyword{internal}
