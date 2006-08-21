## Methods for the sparse QR decomposition

setMethod("qr.qy", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.qy", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.qy", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y) .Call(sparseQR_qty, qr, as.matrix(as.double(y)),
                                FALSE, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.qty", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, TRUE, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.qty", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y) .Call(sparseQR_qty, qr, y, FALSE, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.qty", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y) .Call(sparseQR_qty, qr, as.matrix(as.double(y)),
                                FALSE, TRUE),
          valueClass = "dgeMatrix")

.coef.trunc <- function(qr, res) res[1:ncol(qr@R),,drop=FALSE]

setMethod("qr.coef", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr, y, TRUE)),
          valueClass = "dgeMatrix")

setMethod("qr.coef", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr, y, FALSE)),
          valueClass = "dgeMatrix")

setMethod("qr.coef", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y)
          .coef.trunc(qr, .Call(sparseQR_coef, qr,
                                as.matrix(as.double(y)), FALSE)),
          valueClass = "dgeMatrix")

setMethod("qr.resid", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y)
          .Call(sparseQR_resid_fitted, qr, y, TRUE, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.resid", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y)
          .Call(sparseQR_resid_fitted, qr, y, FALSE, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.resid", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y) .Call(sparseQR_resid_fitted, qr,
                                as.matrix(as.double(y)), FALSE, TRUE),
          valueClass = "dgeMatrix")

setMethod("qr.fitted", signature(qr = "sparseQR", y = "dgeMatrix"),
          function(qr, y, k)
          .Call(sparseQR_resid_fitted, qr, y, TRUE, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.fitted", signature(qr = "sparseQR", y = "matrix"),
          function(qr, y, k)
          .Call(sparseQR_resid_fitted, qr, y, FALSE, FALSE),
          valueClass = "dgeMatrix")

setMethod("qr.fitted", signature(qr = "sparseQR", y = "numeric"),
          function(qr, y, k) .Call(sparseQR_resid_fitted, qr,
                                   as.matrix(as.double(y)), FALSE, FALSE),
          valueClass = "dgeMatrix")
