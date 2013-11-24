\name{rankMatrix}
\Rdversion{1.1}
\alias{rankMatrix}
\title{Rank of a Matrix}
\description{
  Compute \sQuote{the} matrix rank, a well-defined functional in theory,
  somewhat ambigous in practice.  We provide several methods, the
  default corresponding to Matlab's definition.
}
\usage{
rankMatrix(x, tol = NULL,
           method = c("tolNorm2", "qr", "useGrad", "maybeGrad"),
           sval = svd(x, 0, 0)$d, warn.t = TRUE)
}
\arguments{
  \item{x}{numeric matrix, of dimension \eqn{n \times m}{n x m}, say.}
  \item{tol}{nonnegative number specifying a tolerance for
    \dQuote{practically zero} with specific meaning depending on
    \code{method}; by default, \code{max(dim(x)) *
      \link{.Machine}$double.eps * abs(max(sval))} is according to
    Matlab's default (for its only \code{method} "tolNorm2").}
  \item{method}{a character string specifying the computational method
    for the rank, can be abbreviated:
    \describe{
      \item{\code{"tolNorm2"}:}{the number of singular values \code{>= tol};}
      \item{\code{"qr"}:}{for a dense matrix, this is the rank of
	\code{\link[base]{qr}(x, tol, LAPACK=FALSE)} (which is
	\code{qr(...)$rank});

	for sparse \code{x}, where \code{qr()} uses a "sparseQR" method (see
	\code{\link{qr-methods}}) to compute \eqn{QR}, it is the rank of
	\eqn{R} defined as the number of non-"zero" diagonal entries
	\eqn{d_i} of \eqn{R}, and non-"zero"s fulfill \eqn{d_i > tol * max(d_i)}.

	This (dense version) used to be \emph{the} recommended way to
	compute a matrix rank for a while in the past.  For this method,
	\code{sval} are not used (nor computed), which may be crucially
	important for a large sparse matrix \code{x}.}
      \item{\code{"useGrad"}:}{considering the \dQuote{gradient} of the
	(decreasing) singular values, the index of the \emph{smallest} gap.}
      \item{\code{"maybeGrad"}:}{choosing method \code{"useGrad"} only when
	that seems \emph{reasonable}; otherwise using \code{"tolNorm2"}.}

%% FIXME say more
    }
  }
  \item{sval}{numeric vector of non-increasing singular values of
    \code{x}; typically unspecified and computed from \code{x} when
    needed, i.e., unless \code{method = "qr"}.}
  \item{warn.t}{logical indicating if \code{rankMatrix()} should warn
    when it needs \code{\link{t}(x)} instead of \code{x}. Currently, for
    \code{method = "qr"} only.}
}
% \details{
%   FIXME
% }
\note{For large sparse matrices \code{x}, unless you can specify
  \code{sval} yourself, currently \code{method = "qr"} may
  be the only feasible one, as the others need \code{sval} and call
  \code{\link{svd}()} which currently coerces \code{x} to a
  \code{\linkS4class{denseMatrix}} which may be very slow or impossible,
  depending on the matrix dimensions.

  Note that in the case of sparse \code{x}, \code{method = "qr"}, all
  non-strictly zero diagonal entries \eqn{d_i} where counted up to
  including \pkg{Matrix} version 1.1-0, i.e., that method implicitly
  used \code{tol = 0}, see also the seed(42) example below.
}
\value{
  positive integer in \code{1:min(dim(x))}, with attributes detailing
  the method used.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{Martin Maechler; for the "*Grad" methods, building on
  suggestions by Ravi Varadhan.
}
\seealso{
 \code{\link{qr}}, \code{\link{svd}}.
}
\examples{
rankMatrix(cbind(1, 0, 1:3)) # 2

(meths <- eval(formals(rankMatrix)$method))

## a "border" case:
H12 <- Hilbert(12)
rankMatrix(H12, tol = 1e-20) # 12;  but  11  with default method & tol.
sapply(meths, function(.m.) rankMatrix(H12, method = .m.))
## tolNorm2 qr   useGrad maybeGrad
##       11        12        11        11

## "sparse" case:
M15 <- kronecker(diag(x=c(100,1,10)), Hilbert(5))
sapply(meths, function(.m.) rankMatrix(M15, method = .m.))
#--> all 15, but 'useGrad' has 14.

## "large" sparse
n <- 200000; p <- 33; nnz <- 10000
L <- sparseMatrix(i = sample.int(n, nnz, replace=TRUE),
                  j = sample.int(p, nnz, replace=TRUE), x = rnorm(nnz))
(st1 <- system.time(r1 <- rankMatrix(L)))                # almost 1 sec
(st2 <- system.time(r2 <- rankMatrix(L, method = "qr"))) # considerably faster!
r1[[1]] == print(r2[[1]]) ## -->  ( 33  TRUE )

## another sparse-"qr" one, which ``failed'' till 2013-11-23:
set.seed(42)
f1 <- factor(sample(50, 1000, replace=TRUE))
f2 <- factor(sample(50, 1000, replace=TRUE))
f3 <- factor(sample(50, 1000, replace=TRUE))
D <- t(do.call('rBind', lapply(list(f1,f2,f3), as, 'sparseMatrix')))
dim(D); nnzero(D) ## 1000 x 150 // 3000 non-zeros (= 2\%)
stopifnot(rankMatrix(D,           method='qr') == 148,
	  rankMatrix(crossprod(D),method='qr') == 148)
}
\keyword{algebra}
\keyword{array}
