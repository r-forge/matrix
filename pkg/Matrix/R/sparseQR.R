## METHODS FOR CLASS: sparseQR
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## TODO: define implicit generic qr.(Q|R|X|qy|qty|coef|fitted|resid) with
##       formal argument '...' so that we can define methods with further
##       optional arguments ('backPermute', etc.) and perhaps deprecate
##       the 'qrR' work-around

setMethod("expand2", "sparseQR",
          function(x, complete = FALSE, ...) {
              d <- x@Dim
              m <- d[1L]
              n <- d[2L]
              r <- if(complete) m else n
              dn <- x@Dimnames
              p1 <- x@p
              p2 <- x@q
              P1. <- new("pMatrix",
                         Dim = c(m, m),
                         Dimnames = c(dn[1L], list(NULL)),
                         margin = 1L,
                         perm = invPerm(p1, zero.p = TRUE, zero.res = FALSE))
              P2. <- new("pMatrix",
                         Dim = c(n, n),
                         Dimnames = c(list(NULL), dn[2L]),
                         margin = 2L,
                         perm = if(length(p2)) invPerm(p2, zero.p = TRUE, zero.res = FALSE) else seq_len(n))
              ## FIXME: C code should not copy 'y' if it is unreferenced
              ##        and should have an option to _not_ permute !!
              Q <- .Call(sparseQR_qty, x, diag(x = 1, nrow = m, ncol = r),
                         FALSE, FALSE)[p1 + 1L, , drop = FALSE]
              R <- x@R
              if(R@Dim[1L] > r)
                  R <- R[seq_len(r), , drop = FALSE]
              if(complete)
                  list(P1. = P1., Q = Q, R = R, P2. = P2.)
              else
                  list(P1. = P1., Q1 = Q, R1 = triu(R), P2. = P2.)
          })

qrR <- function(qr, complete = FALSE, backPermute = TRUE, row.names = TRUE) {
    d <- qr@Dim
    m <- d[1L]
    n <- d[2L]
    r <- if(complete) m else n
    dn <- qr@Dimnames
    p2 <- qr@q + 1L
    p2.uns <- is.unsorted(p2, strictly = TRUE) # FALSE if length is 0
    R <- qr@R
    m0 <- R@Dim[1L]
    if(!row.names)
        dn <- c(list(NULL), dn[2L])
    else if(m0 > m && !is.null(dn[[1L]]))
        length(dn[[1L]]) <- m0
    if(p2.uns && !is.null(cn <- dn[[2L]]))
        dn[[2L]] <- cn[p2]
    R@Dimnames <- dn
    R <-
        if(m0 > r) {
            if(backPermute && p2.uns)
                R[seq_len(r), invPerm(p2), drop = FALSE]
            else R[seq_len(r), , drop = FALSE]
        } else {
            if(backPermute && p2.uns)
                R[, invPerm(p2), drop = FALSE]
            else R
        }
    if(complete || backPermute) R else triu(R)
}

setMethod("qr.Q", "sparseQR",
	  function(qr, complete = FALSE, Dvec) {
              d <- qr@Dim
              m <- d[1L]
              r <- if(complete) m else d[2L]
              .Call(sparseQR_qty, qr, diag(x = 1, nrow = m, ncol = r),
                    FALSE, FALSE)
          })

setMethod("qr.R", signature(qr = "sparseQR"),
	  function(qr, complete = FALSE) {
              if(nonTRUEoption("Matrix.quiet.qr.R") &&
                 nonTRUEoption("Matrix.quiet"))
		  warning("qr.R(qr(<sparse>)) may differ from qr.R(qr(<dense>)) due to pivoting; see help(\"qr-methods\") and help(\"sparseQR-class\")")
	      qrR(qr, complete = complete, backPermute = FALSE)
	  })

setMethod("qr.qy", signature(qr = "sparseQR", y = "ddenseMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE, TRUE))

setMethod("qr.qy", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE, TRUE))

setMethod("qr.qy", signature(qr = "sparseQR", y = "numeric"),
	  function(qr, y) .Call(sparseQR_qty, qr, y, FALSE, TRUE)@x)

setMethod("qr.qy", signature(qr = "sparseQR", y = "Matrix"),
	  function(qr, y) .Call(sparseQR_qty, qr, .dense2g(as(y, "denseMatrix"), "d"), FALSE, TRUE))

setMethod("qr.qty", signature(qr = "sparseQR", y = "ddenseMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE, TRUE))

setMethod("qr.qty", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE, TRUE))

setMethod("qr.qty", signature(qr = "sparseQR", y = "numeric"),
	  function(qr, y) .Call(sparseQR_qty, qr, y, TRUE, TRUE)@x)

setMethod("qr.qty", signature(qr = "sparseQR", y = "Matrix"),
	  function(qr, y) .Call(sparseQR_qty, qr, .dense2g(as(y, "denseMatrix"), "d"), TRUE, TRUE))

## FIXME: really should happen in C {in sparseQR_coef() in ../src/sparseQR.c} :
.coef.trunc <- function(qr, res) {
    ldn <- lengths(res@Dimnames, use.names = FALSE)
    d <- res@Dim
    m <- qr@R@Dim[2L]
    if(any(ldn > 0L | ldn != d)) {
	## Fix dimnames from dim (when not NULL !) :
	if(ldn[1L]) length(res@Dimnames[[1L]]) <- d[1L]
	if(ldn[2L]) length(res@Dimnames[[2L]]) <- d[2L]
    }
    if(d[1L] == m) res else res[seq_len(m), , drop = FALSE]
}

setMethod("qr.coef", signature(qr = "sparseQR", y = "ddenseMatrix"),
          function(qr, y)
              .coef.trunc(qr, .Call(sparseQR_coef, qr, y)))

setMethod("qr.coef", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              .coef.trunc(qr, .Call(sparseQR_coef, qr, y)))

setMethod("qr.coef", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y)
         drop(.coef.trunc(qr, .Call(sparseQR_coef, qr, y))))

setMethod("qr.coef", signature(qr = "sparseQR", y = "Matrix"),
	  function(qr, y)
              .coef.trunc(qr, .Call(sparseQR_coef, qr, .dense2g(as(y, "denseMatrix"), "d"))))

setMethod("qr.fitted", signature(qr = "sparseQR", y = "ddenseMatrix"),
          function(qr, y, k = qr$rank)
              .Call(sparseQR_resid_fitted, qr, y, FALSE))

setMethod("qr.fitted", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y, k = qr$rank)
              .Call(sparseQR_resid_fitted, qr, y, FALSE))

setMethod("qr.fitted", signature(qr = "sparseQR", y = "numeric"),
	  function(qr, y, k = qr$rank)
         drop(.Call(sparseQR_resid_fitted, qr, y, FALSE)))

setMethod("qr.fitted", signature(qr = "sparseQR", y = "Matrix"),
	  function(qr, y, k = qr$rank)
	      .Call(sparseQR_resid_fitted, qr, .dense2g(as(y, "denseMatrix"), "d"), FALSE))

setMethod("qr.resid", signature(qr = "sparseQR", y = "ddenseMatrix"),
          function(qr, y)
              .Call(sparseQR_resid_fitted, qr, y, TRUE))

setMethod("qr.resid", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
              .Call(sparseQR_resid_fitted, qr, y, TRUE))

setMethod("qr.resid", signature(qr = "sparseQR", y = "numeric"),
	  function(qr, y)
         drop(.Call(sparseQR_resid_fitted, qr, y, TRUE)))

setMethod("qr.resid", signature(qr = "sparseQR", y = "Matrix"),
	  function(qr, y)
              .Call(sparseQR_resid_fitted, qr, .dense2g(as(y, "denseMatrix"), "d"), TRUE))
