\name{mkRespMod}
\alias{mkRespMod}
\title{Create a respModule object}
\description{
  Create a \code{\linkS4class{respModule}} object, which could be from a
  derived class such as \code{\linkS4class{glmRespMod}} or
  \code{\linkS4class{nlsRespMod}}.
}
\usage{
mkRespMod(fr, family = NULL, nlenv = NULL, nlmod = NULL)
}
\arguments{
  \item{fr}{
    a model frame, usually created by a call to \code{\link{model.frame}}.}
  \item{family}{
    an optional glm \code{\link{family}} object
    (\code{\linkS4class{glmRespMod}} objects only).}
  \item{nlenv}{
    an environment for evaluation of the nonlinear model, \code{nlmod}.
    (\code{\linkS4class{nlsRespMod}} objects only).}
  \item{nlmod}{
    the nonlinear model function, as a function call
    (\code{\linkS4class{nlsRespMod}} objects only).}
}
\details{
  The internal representation of a statistical model based on
  a linear predictor expression is derived from a \code{\link{formula}}
  expression and a \code{data} argument, possibly supplemented with a
  \code{\link{family}} object and/or a nonlinear model expression.  The
  steps to obtain this representation usually involve calls to
  \code{\link{model.frame}} and to \code{\link{model.matrix}} or
  \code{\link{model.Matrix}}, which encapsulate important parts of this
  process.   This function encapsulates other operations related to
  weights and offsets and to the model family to create a
  \code{\linkS4class{respModule}} object.
}
\value{
  an object of a class inheriting from \code{\linkS4class{respModule}}.
}
%\references{}
%\author{}
%\note{}
\seealso{
  The \code{\linkS4class{respModule}} class description.
}
\examples{
data(CO2, package="datasets")
## dense linear model
str(Matrix:::glm1(uptake ~ 0 + Type*Treatment, data=CO2, doFit = FALSE))
## sparse linear model
str(Matrix:::glm1(uptake ~ 0 + Type*Treatment, data=CO2, sparse = TRUE, doFit = FALSE))
## From example(glm)
str(trial <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                        outcome=gl(3,1,9,labels=LETTERS[1:3]),
                        treatment=gl(3,3,labels=letters[1:3])))
## dense generalized linear model
str(Matrix:::glm1(counts ~ outcome + treatment, poisson, trial, verbose = TRUE))
## sparse generalized linear model
str(Matrix:::glm1(counts ~ outcome + treatment, poisson, trial, verbose
= TRUE, sparse = TRUE))
# A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                       lot1 = c(118,58,42,35,27,25,21,19,18),
                       lot2 = c(69,35,26,21,18,16,13,12,12))
str(Matrix:::glm1(lot1 ~ log(u), data=clotting, family=Gamma, verbose=TRUE))
}
\keyword{models}
