

## We will no longer use a new 'groupedData' class.  Instead, the
## groupedData function will just add some attributes to the
## data.frame.  The generic function gplot (defined in lattice) will
## be used to produce Trellis plots.






## new version of groupedData

groupedData <-
    function(formula, data,
             order.groups = TRUE,
             FUN = function(x, ...) max(x, na.rm = TRUE),
             outer = NULL, inner = NULL,
             labels = list(), units = list(), ...)

### outer: only used to order panels - same level of outer clustered
### inner: equivalent of groups. Not sure what more than one means
### ...: stored separately, can be used to override gplotArgs calculations
{
    attr(data, "ginfo") <-
        list(formula = formula,
             order.groups = order.groups,
             FUN = FUN,
             outer = outer,
             inner = inner,
             labels = labels,
             units = units)
    dots <- list(...)
    if (length(dots) > 0)
        attr(data, "gplot.args") <- dots
    data
}




convertFromS3groupedData <- function(from)
{
    if (!inherits(from, "groupedData")) return(from)

    data <- from
    for (nm in names(data))
        if (is.factor(data[[nm]]))
            data[[nm]] <- as.character(data[[nm]])
    data <- do.call("data.frame", data)

    formula <- attr(from, "formula")
    old.labels <- attr(from, "labels")
    old.units <- attr(from, "units")

    yvar <- lattice:::.responseName(formula)
    xvar <- lattice:::.covariateName(formula)

    labels <- list()
    units <- list()

    if ("x" %in% names(old.labels)) labels[[xvar]] <- old.labels$x
    if ("y" %in% names(old.labels)) labels[[yvar]] <- old.labels$y
    if ("x" %in% names(old.units)) units[[xvar]] <- old.units$x
    if ("y" %in% names(old.units)) units[[yvar]] <- old.units$y

    order.groups <- attr(from, "order.groups")
    if (is.null(order.groups)) order.groups <- TRUE

    FUN <- attr(from, "FUN")
    if (is.null(FUN)) FUN <- function(x, ...) max(x, na.rm = TRUE)

    groupedData(formula,
                data = data,
                order.groups = order.groups,
                FUN = FUN, 
                outer = attr(from, "outer"),
                inner = attr(from, "inner"),
                labels = labels,
                units = units)
}












