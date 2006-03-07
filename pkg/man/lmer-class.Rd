\name{lmer-class}
\title{Mixed Model Representations and *mer Methods}
\docType{class}
\alias{lmer-class}
\alias{mer-class}
\alias{summary.mer-class}
\alias{summary.lmer-class}
\alias{lmer.ranef-class}
% Methods
\alias{anova,mer-method}
\alias{coef,mer-method}
\alias{coerce,lmer,dsTMatrix-method}
\alias{coerce,mer,dtCMatrix-method}
\alias{coerce,mer,pMatrix-method}
\alias{chol,lmer,ANY-method}
\alias{confint,mer-method}
\alias{deviance,mer-method}
\alias{fitted,mer-method}
\alias{formula,mer-method}
\alias{logLik,mer-method}
%\alias{plot,lmer.coef-method}
\alias{plot,lmer.ranef-method}
\alias{resid,lmer-method}
\alias{resid,mer-method}
\alias{residuals,lmer-method}
\alias{residuals,mer-method}
\alias{show,lmer-method}
\alias{show,mer-method}
\alias{simulate,mer-method}
\alias{summary,mer-method}
\alias{terms,lmer-method}
\alias{update,mer-method}
\alias{vcov,mer-method}
\alias{with,lmer-method}
%
\alias{deviance,summary.mer-method}
\alias{logLik,summary.mer-method}
\alias{vcov,summary.mer-method}
\alias{summary,summary.mer-method}
\description{
  The \code{mer} class is mixed-effects representation using a
  blocked, sparse, symmetric, column-oriented matrix and dense
  matrices.  The \code{lmer} class extends \code{mer} to represent a
  fitted linear mixed-effects model.  The \code{summary.mer} and
  \code{summary.lmer} classes represent summaries of these objects.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lmer", ...)}
  or, more commonly via the \code{\link{lmer}} function.%
  %% Doug, will we want to keep this?
  %% or by applying the \code{summary} generic to an \code{lmer} object.
}

\section{Slots}{
  \describe{
    \item{\code{flist}:}{The list of (possibly permuted) grouping
      factors for the random effects.}
    \item{\code{perm}:}{A list of permutation vectors applied to the
      grouping factors.}
    \item{\code{Parent}:}{The extended parent array for the diagonal
      elements of \code{L}.}
    \item{\code{bVar}:}{A list of the diagonal inner blocks (upper
      triangles only) of the positive-definite matrices on the diagonal
      of the inverse of Z'Z+W.}
    \item{\code{L}:}{The blocked sparse representation of the unit lower
      triangular matrix in the LDL' factorization of Z'Z+W.}
    \item{\code{ZZpO}:}{Object of class \code{"list"} ~~ }
    \item{\code{Omega}:}{A list of numeric matrices providing the
      components of symmetric, positive-definite matrix Omega.  Only the
      upper triangle of each component is used and stored.}
    \item{\code{D}:}{A list of the diagonal factors (upper triangle) in the
      LDL' decomposition of Z'Z+W.}
    \item{\code{RXX}:}{A matrix which is the (augmented) RXX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{RZX}:}{A matrix which is the (augmented) RZX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{ZtZ}:}{The blocks in the blocked sparse symmetric
      representation of the original Z'Z matrix}
    \item{\code{XtX}:}{The original X'X matrix}
    \item{\code{ZtX}:}{The original Z'X matrix}
    \item{\code{cnames}:}{Column names of the model matrices.}
    \item{\code{devComp}:}{Numeric vector of length 4 giving the
      components used to calculate the deviance. This slot's contents
      are current if \code{status$factored} is \code{TRUE}.}
    \item{\code{deviance}:}{Numeric vector of length 2 containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria.  This slot's contents are current if
      \code{status$factored} is \code{TRUE}.}
    \item{\code{nc}:}{Integer vector of the number of columns in
      (augmented) model matrices and the number of observations in the
      model frame.}
    \item{\code{Gp}:}{Integer vector of 0-based row numbers of the
      \code{ZtX} and \code{RZX} matrices where new groups begin.}
    \item{\code{status}:}{Logical vector of length 2 indicating if the
      object has been factored and if the factorization has been
      inverted.}
    \item{\code{call}:}{A copy of the matched call to the function that
      created the object.}
    \item{\code{useScale}:}{logical indicator of the use of the scale
      factor when calculating variance-covariance objects.  Only present
      in the \code{summary.lmer} class.}
    \item{\code{showCorrelation}:}{Logical indicator of whether to
      display the correlation of the fixed-effects estimates when
      showing the object. Only present in the \code{summary.lmer} class.}
  }
  The \code{"summary.mer"} and \code{"summary.lmer"} classes \emph{contain}
  the \code{"mer"} and \code{"lmer"} classes respectively and have
  additional slots,
  \describe{
    \item{\code{isG}:}{logical indicating if we have a
      \bold{g}eneralized linear mixed effects model.}
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
\S4method{logLik}{mer}(object, REML = object@method == "REML", \dots)
\S4method{vcov}{mer}(object, REML = object@method == "REML",
       useScale = object@useScale, \dots)
}
\arguments{
  \item{object}{object of class \code{mer} or \code{lmer}.}
  \item{REML}{logical indicating if REML should be used.}
  \item{useScale}{logical indicating if a scale should be computed and
    used; otherwise \code{sc <- 1} is used.}
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
logLik(fm2)
(V2 <- vcov(fm2))
}
\keyword{classes}
