\name{condest}
\alias{condest}
\alias{onenormest}

\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
condest(A, t = min(n, 5), normA = norm(A, "1"), silent = FALSE, quiet = TRUE)

onenormest(A, t = min(n, 5), A.x, At.x, n,
           silent = FALSE, quiet = silent,
           iter.max = 10, eps = 4 * .Machine$double.eps)

}
\arguments{
  \item{A}{a square matrix, optional for \code{onenormest()}, where
    instead of \code{A}, \code{A.x} and \code{At.x} can be specified,
    see there.}
  \item{t}{number of columns to use in the iterations.}
  \item{normA}{number; (an estimate of) the 1-norm of \code{A}, by
    default \code{\link{norm}(A, "1")}; may be replaced by an estimate.}
  \item{silent}{logical indicating if warning and (by default)
    convergence messages should be displayed.}
  \item{quiet}{logical indicating if convergence messages should be
    displayed.}
  \item{A.x, At.x}{when \code{A} is missing, these two must be given as
    functions which compute \code{A \%\% x}, or \code{t(A) \%\% x},
    respectively.}
  \item{n}{\code{ == nrow(A)}, only needed when \code{A} is not specified.}
  \item{iter.max}{maximal number of iterations for the 1-norm estimator.}
  \item{eps}{the relaive change that is deemed irrelevant.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{Both functions return a \code{\link{list}};
  \code{onenormest()} with components,
  \item{est}{a number \eqn{> 0}, the estimated \code{norm(A, "1")}.}
  \item{v}{the maximal \eqn{A X} column.}

  The function \code{condest()} returns a list with components,
  \item{est}{a number \eqn{> 0}, the estimated condition number
    \eqn{\hat\kappa}{k.}; when \eqn{r :=}\code{rcond(A)},
    \eqn{1/\hat\kappa \approx r}{1/k. ~= r}.}
  \item{v}{integer vector length \code{n}, with an \code{1} at the index
    \code{j} with maximal column \code{A[,j]} in \eqn{A}.}
  \item{w}{numeric vector, the largest \eqn{A x} found.}
  \item{iter}{the number of iterations used.}
}
\references{

}
%\author{Martin Maechler}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
