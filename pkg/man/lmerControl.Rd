% $Id: lmeControl.Rd,v 1.5 2003/09/20 11:54:43 saikat Exp $
\name{lmerControl}
\alias{lmerControl}
\title{Control values for lme}
\description{
  The values supplied in the function call replace the defaults and a
  list with all possible arguments is returned. The returned list is
  used as the \code{control} argument in the \code{lme} function.
}
\synopsis{
lmerControl(maxIter = 50,
            msMaxIter = 50,
            niterEM = 20,
            msTol = sqrt(.Machine$double.eps),
            msVerbose = getOption("verbose"),
            PQLmaxIt = 20,
            EMverbose = getOption("verbose"),
            analyticGradient = TRUE,
            analyticHessian=FALSE))
}
\usage{
lmerControl(maxIter, msMaxIter, niterEM, msTol,
           msVerbose, PQLmaxIt, EMverbose,
           analyticGradient, analyticHessian)
}
\arguments{
}
\value{
  a list with a component for each of the possible arguments.
}
\seealso{\code{\link{lmer}}, \code{\link[stats]{optim}}}
\examples{
# decrease the maximum number iterations in the ms call and
# request that information on the evolution of the ms iterations be printed
str(lmerControl(msMaxIter = 20, msVerbose = TRUE))
}
\keyword{models}
