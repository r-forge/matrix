## METHODS FOR GENERIC: sort
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.V.sort <-
function(x, decreasing = FALSE, na.last = NA, ...) {
    nnz <- length(x@i)
    if (nnz == 0L)
        return(x)
    n <- length(x)
    kind <- .M.kind(x)
    if (kind == "n") {
        x@i <- if (decreasing)
                   seq_len(nnz)
               else seq.int(to = n, length.out = nnz)
        return(x)
    }
    x@x <- y <- sort.int(x@x, na.last = na.last,
                         decreasing = decreasing, ...)
    if (!is.na(na.last))
        nna <- if (anyNA(y)) sum(is.na(y)) else 0L
    else {
        x@length <- n <- n - (nnz - length(y))
        nna <- 0L
        nnz <- length(y)
    }
    nnn <- switch(kind,
                  "l" = nnz - nna,
                  "i" = sum(y >= 0L, na.rm = TRUE),
                  "d" = sum(y >= 0 , na.rm = TRUE),
                  "z" =
                      {
                          arg <- Arg(y)
                          hpi <- 0.5 * pi
                          sum(arg > -hpi & arg <= hpi, na.rm = TRUE)
                      },
                  stop("should never happen ..."))
    if (nna > 0L && decreasing != na.last)
        nnn <- nnn + nna
    x@i <-
        if (nnn < nnz) {
            if (decreasing)
                c(seq_len(nnn), seq.int(to = n, length.out = nnz - nnn))
            else
                c(seq_len(nnz - nnn), seq.int(to = n, length.out = nnn))
        } else {
            if (decreasing)
                seq_len(nnn)
            else
                seq.int(to = n, length.out = nnn)
        }
    x
}

setMethod("sort", c(x = "sparseVector"), .V.sort)
