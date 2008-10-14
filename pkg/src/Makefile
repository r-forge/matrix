#
pkg = Matrix

PKG_CFLAGS = -I./UFconfig
## we use the BLAS and the LAPACK library:
PKG_LIBS = $(SUBLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

ifdef build_under_Windows
  ALL_LIBS = $(PKG_LIBS) -lR
  PKG_HOME = $(shell $(RHOME)/src/gnuwin32/Rpwd.exe .)
  MkInclude = $(PKG_HOME)/Win.mk
else
  MkInclude = $(R_HOME)/etc${R_ARCH}/Makeconf
endif

include SOURCES_C.mkf

OBJECTS = $(SOURCES_C:.c=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = SPQR CHOLMOD COLAMD AMD
SUBLIBS = $(SUBDIRS:=.a)

## Link with C++
SHLIB_LD  = $(SHLIB_CXXLD)

$(SHLIB): sublibs $(OBJECTS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

sublibs:
	@for d in $(SUBDIRS); do \
	  (cd $${d} && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" MkInclude="$(MkInclude)" $(MAKE) library) || exit 1; \
	done

clean:
	@-rm -rf .libs _libs
	@-rm -f *.o $(SHLIB) *.a
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude="$(MkInclude)" $(MAKE) clean) || exit 1; \
	done

include DEPS.mkf
