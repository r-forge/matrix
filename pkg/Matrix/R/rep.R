## METHODS FOR GENERIC: rep
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.V.rep.each <-
function(x, each) {
    each <- as.double(each)
    if (length(each) != 1L) {
        warning(gettextf("first element used of '%s' argument",
                         "each"),
                domain = NA)
        each <- each[1L]
    }
    if (!is.finite(each) || each <= -1)
        stop(gettextf("invalid '%s' argument",
                      "each"),
             domain = NA)
    if (each < 1)
        return(x[0L])
    if (each < 2)
        return(x)
    n <- length(x)
    each <- trunc(each)
    if (n * each > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s",
                      "sparseVector", "2^53"),
             domain = NA)
    if (n * each > .Machine$integer.max) {
        a <- as.double
        one <- 1
    } else {
        each <- as.integer(each)
        a <- as.integer
        one <- 1L
    }
    x@length <- n * each
    x@i <- rep(each * (a(x@i) - one), each = each) + seq_len(each)
    if (.M.kind(x) != "n")
        x@x <- rep(x@x, each = each)
    x
}

.V.rep.int <-
function(x, times) {
    times <- as.double(times)
    if (length(times) != 1L) {
        ## FIXME: support length(times) == length(x)
        warning(gettextf("first element used of '%s' argument",
                         "times"),
                domain = NA)
        times <- times[1L]
    }
    if (!is.finite(times) || times <= -1)
        stop(gettextf("invalid '%s' argument",
                      "times"),
             domain = NA)
    if (times < 1)
        return(x[0L])
    if (times < 2)
        return(x)
    n <- length(x)
    times <- trunc(times)
    if (n * times > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s",
                      "sparseVector", "2^53"),
             domain = NA)
    if (n * times > .Machine$integer.max) {
        a <- as.double
        zero <- 0
    } else {
        times <- as.integer(times)
        a <- as.integer
        zero <- 0L
    }
    x@length <- n * times
    x@i <- rep(a(seq.int(from = zero, by = n, length.out = times)),
               each = length(x@i)) + x@i
    if (.M.kind(x) != "n")
        x@x <- rep.int(x@x, times)
    x
}

.V.rep.len <-
function(x, length.out) {
    length.out <- as.double(length.out)
    if (length(length.out) != 1L) {
        warning(gettextf("first element used of '%s' argument",
                         "length.out"),
                domain = NA)
        length.out <- length.out[1L]
    }
    if (!is.finite(length.out) || length.out <= -1)
        stop(gettextf("invalid '%s' argument",
                      "length.out"),
             domain = NA)
    if (length.out > 0x1p+53)
        stop(gettextf("%s length cannot exceed %s",
                      "sparseVector", "2^53"),
             domain = NA)
    n <- length(x)
    length.out <-
        if (length.out - 1 < .Machine$integer.max)
            as.integer(length.out)
        else trunc(length.out)
    if (length.out > n && n > 0L) {
        x <- .V.rep.int(x, ceiling(length.out / n))
        n <- length(x)
    }
    x@length <- length.out
    if (length.out < n) {
        head <- x@i <= length.out
        x@i <- x@i[head]
        if (.M.kind(x) != "n")
            x@x <- x@x[head]
    } else if (length.out > n && n == 0L) {
        x@i <- seq_len(length.out)
        if (.M.kind(x) != "n")
            x@x <- rep.int(x@x[NA_integer_], length.out)
    }
    x
}

setMethod("rep", c(x = "denseMatrix"),
          function(x, ...) rep(.M2v(x), ...))

setMethod("rep", c(x = "sparseMatrix"),
          function(x, ...) rep(.M2V(x), ...))

setMethod("rep", c(x = "sparseVector"),
          function(x, times, length.out, each, ...) {
              if (!missing(each))
                  x <- .V.rep.each(x, each)
              if (!missing(length.out))
                  x <- .V.rep.len(x, length.out)
              else if (!missing(times))
                  x <- .V.rep.int(x, times)
              x
          })
