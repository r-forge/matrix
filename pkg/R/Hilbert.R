Hilbert <- function(n)
{   ## generate the Hilbert matrix of dimension n
    n <- as.integer(n)
    i <- seq(length = n)
    new("poMatrix", x = c(1/outer(i - 1, i, "+")), Dim = c(n,n))
}
