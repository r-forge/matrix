setMethod("gsummary", signature(object = "data.frame", groups = "factor"),
          function (object,
                    FUN = function(x) mean(x, na.rm = TRUE),
                    omitGroupingFactor = FALSE, 
                    form = formula(object),
                    level,
                    groups = getGroups(object, form, level),
                    invariantsOnly = FALSE, ...)
      {
          groups <- drop(groups)
          gunique <- levels(groups)
          firstInGroup <- match(gunique, groups)
          asFirst <- firstInGroup[match(groups, gunique)]
          value <- as.data.frame(object[firstInGroup, , drop = FALSE])
          row.names(value) <- as.character(gunique)
          value <- value[as.character(sort(gunique)), , drop = FALSE]
          varying <- unlist(lapply(object, function(column, frst) {
              aux <- as.character(column)
              any(!identical(aux, aux[frst]))
          }, frst = asFirst))
          if (any(varying) && (!invariantsOnly)) {
              Mode <- function(x) names(which.max(table(x)))
              if (is(FUN, "function")) {
                  FUN <- list(numeric = FUN, ordered = Mode, factor = Mode)
              }
              else {
                  if (!(is.list(FUN) && all(sapply(FUN, data.class) == 
                                            "function"))) {
                      stop("FUN can only be a function or a list of functions")
                  }
                  auxFUN <- list(numeric = mean, ordered = Mode, factor = Mode)
                  aux <- names(auxFUN)[is.na(match(names(auxFUN), names(FUN)))]
                  if (length(aux) > 0) 
                      FUN[aux] <- auxFUN[aux]
              }
              for (nm in names(object)[varying]) {
                  dClass <- data.class(object[[nm]])
                  if (dClass == "numeric") {
                      value[[nm]] <- as.vector(tapply(object[[nm]], 
                                                      groups, FUN[["numeric"]], ...))
                  }
                  else {
                      value[[nm]] <- as.vector(tapply(as.character(object[[nm]]), 
                                                      groups, FUN[[dClass]]))
                      if (inherits(object[, nm], "ordered")) {
                          value[[nm]] <- ordered(value[, nm], levels = levels(object[, 
                                                              nm]))[drop = TRUE]
                      }
                      else {
                          value[[nm]] <- factor(value[, nm], levels = levels(object[, 
                                                             nm]))[drop = TRUE]
                      }
                  }
              }
          }
          else {
              value <- value[, !varying, drop = FALSE]
          }
          if (omitGroupingFactor) {
              if (is.null(form)) {
                  stop("Cannot omit grouping factor without \"form\"")
              }
              grpForm <- getGroupsFormula(form, asList = TRUE)
              if (missing(level)) 
                  level <- length(grpForm)
              grpNames <- names(grpForm)[level]
              whichKeep <- is.na(match(names(value), grpNames))
              if (any(whichKeep)) {
                  value <- value[, whichKeep, drop = FALSE]
              }
              else {
                  return(NULL)
              }
          }
          value
      })




        
