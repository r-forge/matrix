(m6 <- 1 + as(diag(0:5), "dgeMatrix"))
rcond(m6)
solve(m6, Matrix(1:6))
solve(m6, matrix(1:6))    # method dispatch incorrect for integer matrices
solve(m6, matrix(c(1,2,3,4,5,6)))
##solve(m6, 1:6)            # protection stack overflow
