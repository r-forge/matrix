stopifnot(require(Matrix), require(methods)) # Matrix classes; new, slot<-

wrld_1deg <-
    local({
	load(system.file(file.path("external", "wrld_1deg_slots.rda"),
                         package = "Matrix"))
	## -> 'L'
	r <- new("dsCMatrix")
	for (n in c("Dim", "i","p","x"))
	    slot(r, n) <- L[[n]]
	r
    })

## The reverse:
##  L <- list()
##  for (n in c("Dim", "i","p","x"))    L[[n]] <- slot(wrld_1deg, n)

