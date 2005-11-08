#
pkg = Matrix

PKG_CFLAGS = -DUSE_CHOLMOD
PKG_LIBS = ${LAPACK_LIBS} ${BLAS_LIBS}
include SOURCES_C.mkf

## zherk.f is part of BLAS, but not of R's BLAS sources
SOURCES_F = zpotf2.f zpotrf.f zherk.f

OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD UMFPACK COLAMD CCOLAMD AMD Metis LDL
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

## Variables that had been defined in UFconfig/UFconfig.mk
#------------------------------------------------------------------------------
# METIS, used by CHOLMOD
#------------------------------------------------------------------------------

# The path is relative to where it is used, in CHOLMOD/Lib, CHOLMOD/MATLAB, etc.
# You may wish to use an absolute path.
METIS_PATH=../../Metis
METIS = ../../Metis/libmetis.a

#------------------------------------------------------------------------------
# UMFPACK configuration:
#------------------------------------------------------------------------------

# Configuration flags for UMFPACK.  See UMFPACK/Source/umf_config.h for details.
#
# -DNBLAS	do not use the BLAS.  UMFPACK will be very slow.
# -DCBLAS	use the C-BLAS interface to the BLAS
# -DLP64	UMFPACK is being compiled in 64-bit mode
# -DLONGBLAS	use "long" integers for the BLAS
# -DNSUNPERF	do not use the Sun Perf. Library (default is use it on Solaris)
# -DNSCSL	do not use the SGI SCSL BLAS (default is to use it on SGI)
# -DNPOSIX	do not use POSIX routines sysconf and times.
# -DGETRUSAGE	use getrusage
# -DNO_TIMER	do not use any timing routines
# -DNRECIPROCAL	do not multiply by the reciprocal
# -DNO_DIVIDE_BY_ZERO	do not divide by zero

UMFPACK_CONFIG = 

# CHOLMOD Library Modules, which appear in libcholmod.a:
# Core		requires: none
# Check		requires: Core
# Cholesky	requires: Core, AMD, COLAMD.  optional: Partition, Supernodal
# MatrixOps	requires: Core
# Modify	requires: Core
# Partition	requires: Core, CCOLAMD, METIS.  optional: Cholesky
# Supernodal	requires: Core, BLAS, LAPACK
#
# CHOLMOD test/demo Modules (all are GNU GPL, do not appear in libcholmod.a):
# Tcov		requires: Core, Check, Cholesky, MatrixOps, Modify, Supernodal
#		optional: Partition
# Valgrind	same as Tcov
# Demo		requires: Core, Check, Cholesky, MatrixOps, Supernodal
#		optional: Partition
#
# Configuration flags:
# -DNCHECK	    do not include the Check module.	   License GNU LGPL
# -DNCHOLESKY	    do not include the Cholesky module.	   License GNU LGPL
# -DNPARTITION	    do not include the Partition module.   License GNU LGPL
# -DNGPL	    do not include any GNU GPL Modules in the CHOLMOD library:
# -DNMATRIXOPS	    do not include the MatrixOps module.   License GNU GPL
# -DNMODIFY	    do not include the Modify module.      License GNU GPL
# -DNSUPERNODAL     do not include the Supernodal module.  License GNU GPL
#
# -DNPRINT	    do not print anything.
# -D'LONGBLAS=long' or -DLONGBLAS='long long' defines the integers used by
#  		    	LAPACK and the BLAS (defaults to 'int')
# -DNSUNPERF	    for Solaris only.  If defined, do not use the Sun
#			Performance Library

CHOLMOD_CONFIG = 

BLAS = ${BLAS_LIBS}
LAPACK = ${LAPACK_LIBS}
CHOLMOD_CONFIG = ${CPICFLAGS}

## end of UFconfig/UFconfig.mk rules

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) lib) || exit 1; \
	  touch $${d}.stamp; \
	done

$(SUBSTAMP):
	touch $(SUBSTAMP)

clean: 
	@-rm -rf .libs _libs
	@-rm -f *.o *.so *.a $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) clean) || exit 1; \
	  rm -f $${d}.stamp; \
	done

include DEPS.mkf
