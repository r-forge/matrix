## Objects that must be defined early in collation because we use them
## _at top level_ in one or more ./*.R

## Replaces body(F) with result of substituting quoted arguments,
## where typically 'F' is a template for a set of "parallel" methods
updateBody <- function(F, ...) {
    body(F) <- do.call(substitute, list(body(F), as.list(sys.call())[-1L]))
    F
}
