library(Matrix)

### Matrix Products including  cross products

## checking;  'show' is for convenience of the developer
assert.EQ.mat <- function(M, m, tol = if(show) 0 else 1e-15, show=FALSE) {
    ## temporary fix for R-2.0.1
    MM <- as(M, "matrix")
    attr(MM, "dimnames") <- NULL
    if(show) all.equal(MM, m, tol = tol)
    else stopifnot(all.equal(MM, m, tol = tol))
}
## The relative error typically returned by all.equal:
relErr <- function(target, current)
    mean(abs(target - current)) / mean(abs(target))

m5 <- 1 + as(diag(-1:4)[-5,], "dgeMatrix")
## named dimnames:
dimnames(m5) <- list(Rows= LETTERS[1:5], paste("C", 1:6, sep=""))
stopifnot(dim(m5) == 5:6,
          class(cm5 <- crossprod(m5)) == "dpoMatrix")
assert.EQ.mat((c.m5 <- t(m5) %*% m5), as(cm5, "matrix"))
## but the 'dimnames' are not the same (and are *both*) wrong -- FIXME

## right and left "numeric" multiplication
stopifnot(dim(crossprod(t(m5))) == c(5,5),
          class(p1 <- m5 %*% c(10, 2:6)) == "dgeMatrix")
assert.EQ.mat(p1, cbind(c(20,30,33,38,54)))
(p2 <- c(10, 2:5) %*% m5)


proc.time() # for ``statistical reasons''
