library(Matrix)

## Matrix Exponential

assert.EQ.mat <- function(M, m, tol = 1e-15)
    stopifnot(all.equal(as(M, "matrix"), m, tol = tol))


m1 <- Matrix(c(1,0,1,1), nc = 2)
e1 <- expm(m1)
assert.EQ.mat(e1, cbind(c(exp(1),0), exp(1)))

m2 <- Matrix(c(-49, -64, 24, 31), nc = 2)
e2 <- expm(m2)
assert.EQ.mat(e2, rbind(c(-0.73575875814475, 0.55181909965810),
			c(-1.47151759908826, 1.10363824071557)),
	      tol = 1e-14)
## The ``surprising identity''      det(exp(A)) == exp( tr(A) )
## or                           log det(exp(A)) == tr(A) :
stopifnot(all.equal(determinant(e2)$modulus, sum(diag(m2))))

m3 <- Matrix(cbind(0,rbind(6*diag(3),0)), nc = 4)
e3 <- expm(m3)
assert.EQ.mat(e3,
	      rbind(c(1,6,18,36),
		    c(0,1, 6,18),
		    c(0,0, 1, 6),
		    c(0,0, 0, 1)))

