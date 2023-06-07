\name{Cholesky-methods}
\title{Methods for Cholesky Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
\concept{Choleski} % alternate English spelling
%
\alias{Cholesky}
\alias{Cholesky-methods}
%
\alias{Cholesky,ddiMatrix-method}
\alias{Cholesky,diagonalMatrix-method}
\alias{Cholesky,dsCMatrix-method}
\alias{Cholesky,dsRMatrix-method}
\alias{Cholesky,dsTMatrix-method}
\alias{Cholesky,dspMatrix-method}
\alias{Cholesky,dsyMatrix-method}
\alias{Cholesky,generalMatrix-method}
\alias{Cholesky,matrix-method}
\alias{Cholesky,symmetricMatrix-method}
\alias{Cholesky,triangularMatrix-method}
%
\alias{.SuiteSparse_version}
%
\description{
  Computes the pivoted Cholesky factorization of an
  \eqn{n \times n}{n-by-n} real, symmetric matrix \eqn{A},
  which has the general form
  \deqn{P_{1} A P_{1}' = L_{1} D L_{1}' \overset{D_{jj} \ge 0}{=} L L'}{P1 * A * P1' = L1 * D * L1' [ = L * L' ]}
  or (equivalently)
  \deqn{A = P_{1}' L_{1} D L_{1}' P_{1} \overset{D_{jj} \ge 0}{=} P_{1}' L L' P_{1}}{A = P1' L1 * D * L1' * P1 [ = P1' * L * L' * P1 ]}
  where
  \eqn{P_{1}}{P1} is a permutation matrix,
  \eqn{L_{1}}{L1} is a unit lower triangular matrix,
  \eqn{D} is a diagonal matrix, and
  \eqn{L = L_{1} \sqrt{D}}{L = L1 * sqrt(D)}.
  The second equalities hold only for positive semidefinite \eqn{A},
  for which the diagonal entries of \eqn{D} are non-negative
  and \eqn{\sqrt{D}}{sqrt(D)} is well-defined.
  
  Methods for \code{\linkS4class{denseMatrix}} are built on
  LAPACK routines \code{dpstrf}, \code{dpotrf}, and \code{dpptrf}.
  The latter two do not permute rows or columns,
  so that \eqn{P_{1}}{P1} is an identity matrix.

  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CHOLMOD routines \code{cholmod_analyze} and \code{cholmod_factorize_p}.
}
\usage{
Cholesky(A, \dots)
\S4method{Cholesky}{dsyMatrix}(A, perm = TRUE, tol = -1, \dots)
\S4method{Cholesky}{dspMatrix}(A, \dots)
\S4method{Cholesky}{dsCMatrix}(A, perm = TRUE, LDL = !super, super = FALSE,
    Imult = 0, \dots)
\S4method{Cholesky}{ddiMatrix}(A, \dots)
\S4method{Cholesky}{generalMatrix}(A, uplo = "U", \dots)
\S4method{Cholesky}{triangularMatrix}(A, uplo = "U", \dots)
\S4method{Cholesky}{matrix}(A, uplo = "U", \dots)
}
\arguments{
  \item{A}{a \link[=is.finite]{finite}, symmetric, positive
    semidefinite matrix or \code{\linkS4class{Matrix}}
    to be factorized.  If \code{A} is square but not symmetric,
    then it will be \emph{treated} as symmetric; see \code{uplo}.
    Methods for sparse \code{A} do not require positive
    semidefiniteness when \code{LDL = TRUE}, but still signal
    an error for \code{A} whose leading principal minors
    (after pivoting) are zero.}
  \item{perm}{a logical indicating if the rows and columns
    of \eqn{A} should be pivoted.  Methods for sparse \code{A}
    employ the approximate minimum degree (AMD) algorithm
    in order to reduce fill-in, i.e., without regard for
    numerical stability.
    Pivoting for sparsity may introduce nonpositive leading
    principal minors, causing the factorization to fail, in
    which case it may be necessary to set \code{perm = FALSE}.}
  \item{tol}{a \link[=is.finite]{finite} numeric tolerance,
    used only if \code{pivot = TRUE}.
    The factorization algorithm stops if the pivot is less than
    or equal to \code{tol}.  Negative \code{tol} is equivalent
    to \code{nrow(A) * .Machine$double.eps * max(diag(A))}.}
  \item{LDL}{a logical indicating if the simplicial factorization
    should be computed as
    \eqn{P_{1}' L_{1} D L_{1}' P_{1}}{P1' * L1 * D * L1' * P1},
    such that the result stores the diagonal entries of \eqn{D}
    and the strictly lower triangular entries of \eqn{L_{1}}.
    The alternative is \eqn{P_{1}' L L' P_{1}}{P1' * L * L' * P1},
    such that the result stores the lower triangular entries
    of \eqn{L}.  This argument is ignored if \code{super = TRUE}
    (or if \code{super = NA} and the supernodal algorithm is
    chosen), as the supernodal code does not yet support the
    \code{LDL = TRUE} variant.}
  \item{super}{a logical indicating if the factorization should
    use the supernodal algorithm.  The alternative is the simplicial
    algorithm.  Setting \code{super = NA} leaves the choice to
    a CHOLMOD-internal heuristic.}
  \item{Imult}{a \link[=is.finite]{finite} number. The matrix
    that is factorized is \code{A + Imult * diag(nrow(A))},
    i.e., \code{A} plus \code{Imult} times the identity matrix.
    This argument is useful for symmetric, indefinite \code{A},
    as \code{Imult > max(rowSums(abs(A)) - diag(abs(A)))} ensures
    that \code{A + Imult * diag(nrow(A))} is diagonally dominant.
    (Symmetric, diagonally dominant matrices are positive definite.)}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{A} should be used
    to compute the factorization.  The default is \code{"U"},
    even for lower triangular \code{A}, to be consistent with
    \code{\link[base]{chol}} from \pkg{base}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual class \code{\linkS4class{CholeskyFactorization}}.
  For a traditional matrix \code{A}, the specific class is
  \code{\linkS4class{Cholesky}}.
  For \code{A} inheriting from
  \code{\linkS4class{unpackedMatrix}},
  \code{\linkS4class{packedMatrix}}, and
  \code{\linkS4class{sparseMatrix}},
  the specific class is
  \code{\linkS4class{Cholesky}},
  \code{\linkS4class{pCholesky}}, and
  \code{\linkS4class{dCHMsimpl}} or \code{\linkS4class{dCHMsuper}},
  respectively.
}
\details{
  Note that the result of a call to \code{Cholesky} inherits
  from \code{\linkS4class{CholeskyFactorization}} but not
  \code{\linkS4class{Matrix}}.  Users who just want a matrix
  should consider using \code{\link{chol}}, whose methods are
  simple wrappers around \code{Cholesky} returning just the
  upper triangular Cholesky factor \eqn{L'},
  typically as a \code{\linkS4class{triangularMatrix}}.
  An alternative is to extract factors from the list result
  of \code{\link{expand2}(x, LDL=)}, where \code{x} is the
  \code{CholeskyFactorization} object.
  Yet another alternative is to coerce \code{x}, with
  \code{as(x, "dtrMatrix")},
  \code{as(x, "dtpMatrix")}, and
  \code{as(x, "CsparseMatrix")} for \code{x} of class
  \code{\linkS4class{Cholesky}},
  \code{\linkS4class{pCholesky}}, and 
  \code{\linkS4class{dCHMsimpl}} or \code{\linkS4class{dCHMsuper}},
  respectively.  In all cases the result represents \eqn{L}.
  Indeed, in cases where \eqn{L} is needed explicitly but not any
  other factor, coercion should be preferred over \code{expand2}.
}
\seealso{
  Classes \code{\linkS4class{Cholesky}}, \code{\linkS4class{pCholesky}},
  \code{\linkS4class{dCHMsimpl}} and \code{\linkS4class{dCHMsuper}}
  and their methods.

  Classes \code{\linkS4class{dpoMatrix}}, \code{\linkS4class{dppMatrix}},
  and \code{\linkS4class{dsCMatrix}}.

  Generic function \code{\link{chol}},
  for obtaining the upper triangular Cholesky factor \eqn{L'} as a
  \code{\linkS4class{Matrix}}.

  Generic function \code{\link{expand2}},
  for constructing a list of matrix factors from the result.

  Generic functions \code{\link{BunchKaufman}}, \code{\link{Schur}},
  \code{\link{lu}}, and \code{\link{qr}},
  for computing other factorizations.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dpstrf.f},
  \url{https://netlib.org/lapack/double/dpotrf.f}, and
  \url{https://netlib.org/lapack/double/dpptrf.f}.

  Chen, Y., Davis, T. A., Hager, W. W., Rajamanickam, S. (2008).
  Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization
  and update/downdate.
  \emph{ACM Transactions on Mathematical Software},
  \emph{35}(3), Article 22, 1-14.
  \doi{10.1145/1391989.1391995}

  Amestoy, P. R., Davis, T. A., & Duff, I. S. (2004).
  Algorithm 837: AMD, an approximate minimum degree ordering algorithm.
  \emph{ACM Transactions on Mathematical Software},
  \emph{17}(4), 886-905.
  \doi{10.1145/1024074.1024081}
}
\examples{
showMethods("Cholesky", inherited = FALSE)

data(KNex)
mtm <- with(KNex, crossprod(mm))
str(mtm@factors) # empty list()
(C1 <- Cholesky(mtm))             # uses show(<MatrixFactorization>)
str(mtm@factors) # 'sPDCholesky' (simpl)
(Cm <- Cholesky(mtm, super = TRUE))
c(C1 = isLDL(C1), Cm = isLDL(Cm)) # TRUE FALSE
str(mtm@factors) # 'sPDCholesky'  *and* 'SPdCholesky'
str(cm1  <- as(C1, "CsparseMatrix"))
str(cmat <- as(Cm, "CsparseMatrix"))# hmm: super is *less* sparse here
cm1[1:20, 1:20]

b <- matrix(c(rep(0, 711), 1), ncol = 1)
## solve(Cm, b) by default solves  Ax = b, where A = Cm'Cm (= mtm)!
## hence, the identical() check *should* work, but fails on some GOTOblas:
x <- solve(Cm, b)
stopifnot(identical(x, solve(Cm, b, system = "A")),
          all.equal(x, solve(mtm, b)))

Cn <- Cholesky(mtm, perm = FALSE)# no permutation -- much worse:
sizes <- c(simple = object.size(C1),
           super  = object.size(Cm),
           noPerm = object.size(Cn))
## simple is 100, super= 137, noPerm= 812 :
noquote(cbind(format(100 * sizes / sizes[1], digits=4)))


## Visualize the sparseness:
dq <- function(ch) paste('"',ch,'"', sep="") ## dQuote(<UTF-8>) gives bad plots
image(mtm, main=paste("crossprod(mm) : Sparse", dq(class(mtm))))
image(cm1, main= paste("as(Cholesky(crossprod(mm)), \"CsparseMatrix\"):",
                        dq(class(cm1))))
\dontshow{% FIXME-- move to ../tests/factorizing.R :
expand(C1) ## to check printing
}

## Smaller example, with same matrix as in  help(chol) :
(mm <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)) # 5 x 5
(opts <- expand.grid(perm = c(TRUE,FALSE), LDL = c(TRUE,FALSE), super = c(FALSE,TRUE)))
rr <- lapply(seq_len(nrow(opts)), function(i)
             do.call(Cholesky, c(list(A = mm), opts[i,])))
nn <- do.call(expand.grid, c(attr(opts, "out.attrs")$dimnames,
              stringsAsFactors=FALSE,KEEP.OUT.ATTRS=FALSE))
names(rr) <- apply(nn, 1, function(r)
                   paste(sub("(=.).*","\\\\1", r), collapse=","))% extra '\\' in Rd
str(rr, max.level=1)

str(re <- lapply(rr, expand), max.level=2) ## each has a 'P' and a 'L' matrix
%% FIXME !! --- "check" them  __unfinished__
R0 <- chol(mm, pivot=FALSE)
R1 <- chol(mm, pivot=TRUE )
stopifnot(all.equal(t(R1), re[[1]]$L),
          all.equal(t(R0), re[[2]]$L),
          identical(as(1:5, "pMatrix"), re[[2]]$P), # no pivoting
TRUE)

## Symmetric matrix with positive *and* negative eigen values, i.e., indefinite
mpn <- new("dsCMatrix", i = c(0L, 0:1, 0:3, 2:5, 3:6),
           p = c(0:1, 3L, 6:7, 10:11, 15L), Dim = c(7L, 7L),
           x = c(1, 6, 38, 10, 60, 103, -4, 6, -32, -247, -2, -16, -128, -2, -67))
ev <- eigen(mpn, only.values = TRUE)$values
round(ev, 3) # 139.8 3.059 0.268 -0.008 -0.140 -3.357 -317.6
plot(ev, type="b"); abline(h=0, lty=3)
Cpn <- Cholesky(mpn)
diag(Cpn) # D matrix has [p]ositive and [n]egative diagonal elements
try(expand(Cpn)) # cannot compute L = L1 sqrt(D)
(mp. <- as(Reduce(`\%*\%`, expand2(Cpn)), "symmetricMatrix"))
## very close to 'mpn' but *does* differ numerically (!)
stopifnot(all.equal(mpn, mp.))

%% --> ../TODO : .diag.dsC() should be renamed, documented, tested,...
\dontshow{
str(dd <- .diag.dsC(mtm))
dc <- .diag.dsC(Chx=C1) # <- directly from the Cholesky
stopifnot(all.equal(dd,dc))
}%dont
# Version of the underlying SuiteSparse library by Tim Davis :
.SuiteSparse_version()
}
