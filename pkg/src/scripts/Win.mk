## Note: This is included from ./Makefile.win , but also
## ----  as 'MkInclude' from several */Makefile s
include $(RHOME)/src/gnuwin32/MkRules
CFLAGS=$(PKG_CFLAGS) -O3
ALL_CFLAGS=$(PKG_CFLAGS) -O3
CXXFLAGS=$(PKG_CXXFLAGS) -O2 
ALL_CXXFLAGS=$(PKG_CXXFLAGS) -O2
