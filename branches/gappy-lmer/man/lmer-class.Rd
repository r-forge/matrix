\name{lmer-class}
\title{Mixed Model Representations and *mer Methods}
\docType{class}
% Classes
\alias{mer-class}
\alias{lmer-class}
\alias{glmer-class}
\alias{nlmer-class}
\alias{coef.lmer-class}
\alias{summary.mer-class}
\alias{summary.lmer-class}
\alias{summary.glmer-class}
\alias{summary.nlmer-class}
\alias{ranef.mer-class}
% Methods
\alias{anova,lmer-method}
\alias{coef,lmer-method}
\alias{coef,mer-method}
\alias{coerce,lmer,dsTMatrix-method}
\alias{coerce,mer,dtCMatrix-method}
\alias{coerce,mer,pMatrix-method}
\alias{chol,lmer,ANY-method}
\alias{confint,lmer-method}
\alias{deviance,mer-method}
\alias{fitted,lmer-method}
\alias{formula,mer-method}
\alias{logLik,mer-method}
\alias{plot,coef.lmer-method}
\alias{plot,ranef.mer-method}
\alias{print,mer-method}
\alias{print,lmer-method}
\alias{print,glmer-method}
\alias{print,nlmer-method}
\alias{qqmath,ranef.mer-method}
\alias{resid,glmer-method}
\alias{resid,lmer-method}
\alias{residuals,glmer-method}
\alias{residuals,lmer-method}
\alias{show,mer-method}
\alias{show,lmer-method}
\alias{show,glmer-method}
\alias{show,nlmer-method}
\alias{simulate,lmer-method}
\alias{summary,lmer-method}
\alias{summary,mer-method}
\alias{terms,mer-method}
\alias{update,mer-method}
\alias{vcov,mer-method}
\alias{vcov,lmer-method}
\alias{with,lmer-method}
%
\alias{deviance,summary.mer-method}
\alias{logLik,summary.mer-method}
\alias{vcov,summary.mer-method}
\alias{summary,summary.mer-method}
\description{
  The \code{mer} class is mixed-effects representation using a sparse
  model matrix and the sparse Cholesky factor of the (updated)
  cross-product of its columns.  The \code{lmer} class extends
  \code{mer} to represent a fitted linear mixed-effects model.  The
  \code{summary.mer} and \code{summary.lmer} classes represent summaries
  of these objects.  The \code{coef.lmer} and \code{ranef.mer} classes
  represent fixed and random effect \dQuote{coefficients}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lmer", ...)}
  or, more commonly via the \code{\link{lmer}} function.%
  %% Doug, will we want to keep this?
  %% or by applying the \code{summary} generic to an \code{lmer} object.
}

\section{Slots}{
  \describe{
    \item{\code{env}:}{An environment (of class \code{"environment"})
      created for the evaluation of the nonlinear model function.}
    \item{\code{model}:}{The nonlinear model function as an object of
      class \code{"call"}.}
    \item{\code{frame}:}{The model frame for the model, as an object of
      class \code{"data.frame"}.}
    \item{\code{pnames}:}{Names of the parameters in the nonlinear model
      (class \code{"character"}).}
    \item{\code{call}:}{The matched call to the function that
      created the object. (class \code{"call"}).} 
    \item{\code{flist}:}{The list of grouping factors for the random
      effects.}
    \item{\code{Xt}:}{Sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
      the fixed effects.}
    \item{\code{Zt}:}{Sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
      the random effects.}
    \item{\code{y}:}{The response vector (class \code{"numeric"}).}
    \item{\code{weights}:}{Numeric weights vector. This may be of length
      zero (0) indicating unit weights.}
    \item{\code{cnames}:}{a list of character vectors of column names
      of the random-effects variance-covariance matrix associated with
      each grouping factor and the fixed-effects model matrix.}
    \item{\code{Gp}:}{integer vector of group pointers within the random
      effects and the rows of the transposed model matrix in the
      \code{Zt} slot.  The elements of \code{Gp} are the 0-based index
      of the first element corresponding to each grouping factor.  Thus
      the first element is always 0.  The last element is the total
      length of the random effects vector (also the total number of rows
      in the matrix in the \code{Zt} slot).}
    \item{\code{dims}:}{A named integer vector of dimensions.}
    \item{\code{ST}:}{A list of S and T factors in the TSST' Cholesky
      factorization of the relative variance matrices of the random
      effects associated with each grouping factor.  The unit lower
      triangular matrix T and the diagonal matrix S are stored as a
      single matrix whose diagonal elements determine S while the
      subdiagonal elements determine the non-trivial elements of T.}
    \item{\code{Vt}:}{Sparse matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"})
      containing the transpose of \eqn{V=ZTS}.}
    \item{\code{L}:}{The sparse lower Cholesky factor of \eqn{V'V+I}
      (class \code{"\link[Matrix:CHMfactor-class]{dCHMfactor}"}).}
    \item{\code{L}:}{The Cholesky decomposition (class
      \code{"\link[Matrix:CHMfactor-class]{dCHMfactor}"}) of \eqn{V'V+I}.}
    \item{\code{mu}:}{A numeric vector of predicted values from the
      model with an attibute \code{"gradient"}, the gradient matrix.}
    \item{\code{Mt}:}{The derivative of \code{mu} with respect to
      \code{uvec} stored as a sparse matrix (class
      \code{"\link[Matrix:CHMfactor-class]{dCHMfactor}"}).}      
    \item{\code{deviance}:}{Named numeric vector containing the deviance
      corresponding to the maximum likelihood (the \code{"ML"} element)
      and \code{"REML"} criteria and various components. The
      \code{"ldL2"} element is twice the logarithm of the determinant of
      the Cholesky factor in the \code{L} slot. The \code{"lpdisc"}
      element is the logarithm of the penalized discrepancy.  The
      \code{"bqd"} component is the value of the random-effects
      quadratic form.}
    \item{\code{fixef}:}{Numeric vector of fixed effects.}
    \item{\code{ranef}:}{Numeric vector of random effects on the
      original scale.}
    \item{\code{uvec}:}{Numeric vector of orthogonal, constant variance
      random effects.}
    \item{\code{X}:}{Dense model matrix for the fixed
      effects. (\code{"lmer"} class only.)}
    \item{\code{XytXy}:}{The crossproduct of the columns of \eqn{[X:y]}.}
    \item{\code{RVXy}:}{A rectangular matrix, the dense solution to
      \eqn{L RVXy = S T'ZtXy}.}
    \item{\code{RXy}:}{The upper triangular Cholesky factor of the
      downdated \code{XytXy} slot.}
  }
  The \code{"summary.mer"}, \code{"summary.lmer"} and
  \code{"summary.glmer"} classes \emph{contain} the \code{"mer"},
  \code{"lmer"} and \code{"glmer"} classes respectively and have
  additional slots,
  \describe{
    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{\link{logLik}(object)}.}
    \item{\code{ngrps}:}{...}% FIXME
    \item{\code{sigma}:}{the scale ...}% FIXME
    \item{\code{coefs}:}{the coefficients matrix...}% FIXME
    \item{\code{vcov}:}{the same as \code{\link{vcov}(object)}.}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{...}% FIXME
  }% summary slots
}
\usage{
\S4method{logLik}{mer}(object, REML = object@status["REML"], \dots)
\S4method{vcov}{mer}(object, REML = object@status["REML"], \dots)
}
\arguments{
  \item{object}{object of class \code{mer} or \code{lmer}.}
  \item{REML}{logical indicating if REML should be used.}
  \item{\dots}{potential further arguments passed to methods.}
}
\section{Methods}{
% TODO:
%  "show": object = "lmer" (inherited from object = "mer")
%  "with": data = "lmer"
  \describe{
    \item{anova}{\code{signature(object = "mer")}: model comparisons for
      \code{mer} objects, see \code{\link{anova}}.}
      %% FIXME: say a bit more on lmer - model comparison ??
    \item{coef}{\code{signature(object = "lmer")}: Extract the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization.}
    \item{coef<-}{\code{signature(object = "lmer", value = "numeric")}: Assign the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization. }
    \item{coerce}{\code{signature(from = "lmer", to = "dtTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the unit lower triangular matrix \eqn{L}
      from the \eqn{LDL} decomposition of
      \eqn{Z^\prime Z+\Omega}{crossprod(Z)+Omega}.}
    \item{coerce}{\code{signature(from = "lmer", to = "dsTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the symmetric matrix
      \eqn{Z^\prime Z+\Omega}{crossprod(Z)+Omega}.}
    \item{deviance}{\code{signature(object = "mer")}: extracts the
      \code{\link{deviance}} of the fitted model.}
%%% some of the following are also in \usage{..} because they have extra args:
    \item{logLik}{\code{signature(object = "mer")}: Extract the
      log-likelihood from of the fitted model, see also \code{\link{logLik}}.}
    \item{update}{\code{signature(object = "mer")}: see
      \code{\link{update}} on how to update fitted models.}
    \item{simulate}{\code{signature(object = "mer")}: see
      \code{\link{simulate}} on how to simulate from fitted models.}
    \item{solve}{\code{signature(a = "lmer", b = "missing")}: Invert
      the decomposed matrices.}
    \item{terms}{\code{signature(x = "lmer")}: Extract
      \code{\link{terms}} component.}
    \item{vcov}{\code{signature(object = "mer")}: Calculate
      variance-covariance matrix of the \emph{fixed} effect terms,
      see also \code{\link{vcov}}.
    }
  }
}
\seealso{
  \code{\link{lmer}()} produces these objects.
  \cr
  \code{\link{VarCorr}} for extracting the variance and
  correlation components of the \emph{random}-effects terms.
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
print(fm2, digits = 10, corr = FALSE) # more precision; no corr.matrix

logLik(fm2)
(V2 <- vcov(fm2))
\dontrun{
str(resid(fm2))

## Simulate 'Reaction' according to the fitted model:
dim(ss <- simulate(fm2, nsim = 200, seed = 101)) ## -> 180 x 200
}
}
\keyword{classes}
