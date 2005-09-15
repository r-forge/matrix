#-------------------------------------------------------------------------------
# CHOLMOD Makefile
#-------------------------------------------------------------------------------

default: all

# Compile the C-callable libraries and the Demo programs.
all:
	( cd Lib ; make )
	#( cd Demo ; make )

# Compile the C-callable libraries only.
library:
	( cd Lib ; make )

lib: library

# Remove all files not in the original distribution
purge:
	#( cd MATLAB ; make purge )
	( cd Tcov ; make purge )
	( cd Lib ; make purge )
	#( cd Valgrind ; make purge )
	#( cd Demo ; make purge )
	#( cd Doc ; make purge )

# Remove all files not in the original distribution, except keep the 
# compiled libraries.
clean:
	#( cd MATLAB ; make clean )
	( cd Tcov ; make clean )
	( cd Lib ; make clean )
	#( cd Valgrind ; make clean )
	#( cd Demo ; make clean )

distclean: purge

ccode: all

# Compile the MATLAB mexFunctions (you can also use cholmod_make.m in MATLAB)
matlab:
	#( cd MATLAB ; make )

# Run the test coverage suite.  Takes about 20 minutes on a 3.2GHz Pentium.
cov:
	( cd Tcov ; make go )

# Run the test coverage suite using Valgrind.  This takes a *** long *** time.
valgrind:
	( cd Valgrind ; make )

# Compile the C-callable libraries and the Demo programs.
demo:
	#( cd Demo ; make )
