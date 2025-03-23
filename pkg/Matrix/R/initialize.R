## METHODS FOR GENERIC: initialize
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("Matrix", "MatrixFactorization"))
setMethod("initialize", c(.Object = .cl),
          function(.Object, ...) {
              ## Allows Dimnames[[i]] to be a vector of type other than
              ## "character" and to be a vector of length zero instead
              ## of NULL
              .Object <- callNextMethod()
              if (...length() && any(...names() == "Dimnames"))
                  .Object@Dimnames <- fixupDN(.Object@Dimnames)
              .Object
          })
rm(.cl)

setMethod("initialize", c(.Object = "sparseVector"),
          function(.Object, i, x, ...) {
              ## Allows 'i' to be unsorted
              if (!missing(i) &&
                  !is.na(uns <- is.unsorted(i, strictly = TRUE)) &&
                  uns) {
                  ## There are no NA, and the order of ties does not
                  ## matter (since ties are invalid), so it is safe
                  ## to use "quick" here
                  method <- if (is.integer(length(i))) "radix" else "quick"
                  if (missing(x))
                      i <- sort.int(i, method = method)
                  else if (length(x) == length(i)) {
                      s <- sort.int(i, method = method, index.return = TRUE)
                      x <- x[s[["ix"]]]
                      i <- s[["x"]]
                  }
              }
              callNextMethod()
          })
