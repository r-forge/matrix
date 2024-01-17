\name{lmer}
\alias{lmer}
\alias{glmer}
\alias{nlmer}
\concept{GLMM}% << so it's found
\concept{NLMM}% << so it's found
\title{Fit Mixed-Effects Models}
\description{
  Fit a linear mixed model or a generalized linear mixed model or a
  nonlinear mixed model.
}
\usage{
lmer(formula, data, family, method, control, start,
     verbose, subset, weights, na.action, offset,
     contrasts, model, x, \dots)
glmer(formula, data, family, method, control, start,
      verbose, subset, weights, na.action, offset,
      contrasts, model, \dots)
nlmer(formula, data, control, start, verbose,
      subset, weights, na.action, contrasts, model, \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.  The vertical bar character \code{"|"} separates an
    expression for a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called.}
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}}.  If \code{family} is
    missing then a linear mixed model is fit; otherwise a generalized
    linear mixed model is fit.}
  \item{method}{a character string.  For a linear mixed model the
    default is \code{"REML"} indicating that the model should be fit by
    maximizing the restricted log-likelihood.  The alternative is
    \code{"ML"} indicating that the log-likelihood should be
    maximized.  (This method is sometimes called "full" maximum
    likelihood.)  For a generalized linear mixed model the criterion
    is always the log-likelihood but this criterion does not have a
    closed form expression and must be approximated.  The default approximation is
    \code{"PQL"} or penalized quasi-likelihood.  Alternatives are
    \code{"Laplace"} or \code{"AGQ"} indicating the Laplacian and
    adaptive Gaussian quadrature approximations respectively.  The
    \code{"PQL"} method is fastest but least accurate.  The
    \code{"Laplace"} method is intermediate in speed and accuracy.
    The \code{"AGQ"} method is the most accurate but can be
    considerably slower than the others.}
  \item{control}{a list of control parameters. See below for details.}
  \item{start}{a list of relative precision matrices for the random
    effects.  This has the same form as the slot \code{"Omega"} in a
    fitted model.  Only the upper triangle of these symmetric matrices
    should be stored.}
  \item{subset, weights, na.action, offset, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
  \item{model}{logical scalar.  If \code{FALSE} the model frame in
    slot \code{frame} is truncated to zero rows.}
  \item{x}{logical scalar.  If \code{FALSE} the model matrix in
    slot \code{X} is truncated to zero rows.}
  \item{verbose}{logical scalar.  If \code{TRUE} verbose output is
    generated during the optimization of the parameter estimates.}
  \item{\dots}{potentially further arguments for methods.  Currently
    none are used.}
}
\details{
  This is a revised version of the \code{\link[nlme]{lme}} function from the
  \pkg{nlme} package.  This version uses a different method of
  specifying random-effects terms and allows for fitting generalized
  linear mixed models as well as linear mixed models.

  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer}.
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) prints
      an error message and terminate if there are any incomplete
      observations.}
    \item{control}{a named list of control parameters for the estimation
      algorithm, specifying only the ones to be changed from their
      default values.  Hence defaults to an empty list.\cr
      Possible control options and their default values are:
      \describe{
% 	\item{\code{maxIter}:}{maximum number of iterations for the \code{lme}
% 	  optimization algorithm.  Default is 50.}
% 	\item{\code{tolerance}:}{relative convergence tolerance for the
% 	  \code{lme} optimization algorithm.  Default is
% 	  \code{sqrt(.Machine$double.eps)}.}
	\item{\code{msVerbose}:}{a logical value passed as the
	  \code{trace} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{getOption("verbose")}.}
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An object of class \code{"\linkS4class{lmer}"}.
  There are many methods applicable to \code{"lmer"} objects, see the
  above help page. %% Hmm.. or \code{\link{showMethods}(class = "lmer")}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{}
%\note{}
\seealso{
  The \code{\linkS4class{lmer}} class, \code{\link[stats]{lm}}
}
\examples{
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)
(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
              Orange, verb = 1, start = c(Asym = 200, xmid = 725, scal = 350)))
(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKe+lKa+lCl|Subject),
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe|Subject) + (lKa|Subject) + (lCl|Subject), Theoph,
              start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKa+lCl|Subject),
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKa|Subject) + (lCl|Subject), 
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
\dontrun{
(m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp))
}
}
\keyword{models}
\keyword{methods}
