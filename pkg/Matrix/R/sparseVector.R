## METHODS FOR CLASS: sparseVector (virtual)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diff", signature(x = "sparseVector"),
          ## Mostly cut and paste of base::diff.default :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) > 1L ||
                  lag < 1L || differences < 1L)
                  stop("'lag' and 'differences' must be integers >= 1")
              if(lag * differences >= length(x))
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences))
                  x <- x[i1] - x[-length(x):-(length(x) - lag + 1L)]
              x
          })

setMethod("mean", signature(x = "sparseVector"),
          function(x, trim = 0, na.rm = FALSE, ...) {
              kind <- .V.kind(x)
              if(kind == "z" && trim > 0)
                  stop("trimmed means are not defined for complex data")
              n <- length(x)
              if(kind != "n" && n > 0L && anyNA(x@x)) {
                  if(!na.rm)
                      return(NA_real_)
                  n <- n - sum(is.na(x@x))
              }
              if(n == 0L)
                  return(if(kind != "z") NaN else NaN * 0i)
              if(kind == "n") {
                  nnz <- length(x@i)
                  if(trim <= 0)
                      return(nnz / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  if(nnz < ntrim)
                      0
                  else if(nnz == ntrim) {
                      if(n - 2 * ntrim > 0)
                           0
                      else 0.5
                  } else {
                      if(n - 2 * ntrim > 0)
                          (nnz - ntrim - max(ntrim - (n - nnz), 0)) /
                              (n - 2 * ntrim)
                      else 1
                  }
              } else {
                  if(trim <= 0)
                      ## FIXME? not computing in long double
                      return(sum(x@x, na.rm = na.rm) / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  x <- sort(x, na.last = NA)[(ntrim + 1):(n - ntrim)]
                  sum(x@x) / x@length
              }
          })

.V.rep.each <- function(x, each) {
    each <- trunc(as.double(each)[1L])
    if(!is.finite(each) || each < 0)
        stop("invalid 'each' argument")
    if(each == 0)
        return(x[0L])
    if(each == 1)
        return(x)
    n <- length(x)
    if(n * each <= .Machine$integer.max) {
        each <- as.integer(each)
        x@length <- n * each
        x@i <- rep(each * (as.integer(x@i) - 1L), each = each) + seq_len(each)
    } else {
        x@length <- n * each
        x@i <- rep(each * (as.double (x@i) - 1 ), each = each) + seq_len(each)
    }
    if(.hasSlot(x, "x"))
        x@x <- rep(x@x, each = each)
    x
}

.V.rep.int  <- function(x, times) {
    times <- trunc(as.double(times)[1L])
    if(!is.finite(times) || times < 0)
        stop("invalid 'times' argument")
    if(times == 0)
        return(x[0L])
    if(times == 1)
        return(x)
    n <- length(x)
    if(n * times <= .Machine$integer.max) {
        times <- as.integer(times)
        x@length <- n * times
        x@i <- rep(seq.int(from = 0L, by = n, length.out = times),
                   each = length(x@i)) + as.integer(x@i)
    } else {
        x@length <- n * times
        x@i <- rep(seq.int(from = 0 , by = n, length.out = times),
                   each = length(x@i)) + as.double (x@i)
    }
    if(.hasSlot(x, "x"))
        x@x <- rep.int(x@x, times)
    x
}

.V.rep.len  <- function(x, length.out) {
    length.out <- trunc(as.double(length.out)[1L])
    if(!is.finite(length.out) || length.out < 0)
        stop("invalid 'length.out' argument")
    n <- length(x)
    if(length.out > n) {
        if(n == 0L) {
            if(length.out <= .Machine$integer.max)
                length.out <- as.integer(length.out)
            x@length <- length.out
            x@i <- seq_len(length.out)
            if(.hasSlot(x, "x"))
                x@x <- rep.int(x@x[NA_integer_], length.out)
            return(x)
        }
        x <- .V.rep.int(x, ceiling(length.out / n))
        n <- length(x)
    }
    if(length.out == n) x else x[seq_len(length.out)]
}

setMethod("rep", signature(x = "sparseVector"),
          function(x, times, length.out, each, ...) {
              if(!missing(each))
                  x <- .V.rep.each(x, each)
              if(!missing(length.out))
                  x <- .V.rep.len (x, length.out)
              else if(!missing(times))
                  x <- .V.rep.int (x, times)
              x
          })

setMethod("sort", signature(x = "sparseVector"),
          function(x, decreasing = FALSE, na.last = NA, ...) {
              nnz <- length(x@i)
              if(nnz == 0L)
                  return(x)
              n <- length(x)
              kind <- .V.kind(x)
              if(kind == "n") {
                  x@i <- if(decreasing)
                             seq_len(nnz)
                         else seq.int(to = n, length.out = nnz)
                  return(x)
              }
              x@x <- y <- sort.int(x@x, na.last = na.last,
                                   decreasing = decreasing, ...)
              if(!is.na(na.last))
                  nna <- if(anyNA(y)) sum(is.na(y)) else 0L
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
              if(nna > 0L && decreasing != na.last)
                  nnn <- nnn + nna
              x@i <-
                  if(nnn < nnz) {
                      if(decreasing)
                      c(seq_len(nnn), seq.int(to = n, length.out = nnz - nnn))
                      else
                      c(seq_len(nnz - nnn), seq.int(to = n, length.out = nnn))
                  } else {
                      if(decreasing)
                        seq_len(nnn)
                      else
                                            seq.int(to = n, length.out = nnn)
                  }
              x
          })

setMethod("t", signature(x = "sparseVector"),
          function(x) .tCRT(.V2C(x)))

setMethod("toeplitz", signature(x = "sparseVector"),
          function(x, symmetric = TRUE, repr = c("C", "R", "T"),
                   giveCsparse, ...) {
              n <- length(x)
              if(n > .Machine$integer.max)
                  stop("dimensions cannot exceed 2^31-1")
              nn <- c(n, n)
              r <- spV2M(x[as.integer(abs(.col(nn) - .row(nn))) + 1L],
                         nrow = n, ncol = n, symmetric = symmetric,
                         check = FALSE)
              repr <- # keep in sync with sparseMatrix
                  if(missing(giveCsparse))
                      match.arg(repr)
                  else if(!missing(repr)) {
                      warning("'giveCsparse' is deprecated; using 'repr' instead")
                      match.arg(repr)
                  } else if(giveCsparse) {
                      "C"
                  } else {
                      warning("'giveCsparse' is deprecated; setting repr=\"T\" for you")
                      "T"
                  }
              switch(repr, "C" = .M2C(r), "R" = .M2R(r), "T" = r)
          })
