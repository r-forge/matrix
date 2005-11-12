#
pkg = Matrix

PKG_CFLAGS = -DUSE_CHOLMOD
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS)
include SOURCES_C.mkf
## zherk.f is part of BLAS, but not of R's BLAS sources
ifeq ($(strip $(BLAS_LIBS)),)
SOURCES_F = zpotf2.f zpotrf.f zherk.f lsame.f
else
SOURCES_F =
endif

OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD UMFPACK COLAMD CCOLAMD AMD Metis LDL
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

## Variables that had been defined in UFconfig/UFconfig.mk
## (but MM thinks are nowhere used)
BLAS = $(BLAS_LIBS)
LAPACK = $(LAPACK_LIBS)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude=$(R_HOME)/etc/Makeconf $(MAKE) lib) || exit 1; \
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
