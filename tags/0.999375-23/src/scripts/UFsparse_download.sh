#!/bin/sh
## Update Libraries from Tim Davis' University of Florida (UF) collection:
#
if [ ! -d ../src ]
then echo 'Must run in Matrix/src/ !' ; exit 1
fi
ufl_URL=http://www.cise.ufl.edu/research/sparse
wget $ufl_URL/amd/current/AMD.tar.gz
wget $ufl_URL/cholmod/current/CHOLMOD.tar.gz
wget $ufl_URL/colamd/current/COLAMD.tar.gz
wget $ufl_URL/UFconfig/current/UFconfig.tar.gz
wget $ufl_URL/SPQR/current/SPQR.tar.gz

## 1) UFconfig ---------------------------------------------
  ## install UFconfig.h file (now needed by some UFsparse libraries)
tar zxf UFconfig.tar.gz UFconfig/UFconfig.h UFconfig/README.txt
  ## Move the UFconfig/README.txt file to ../inst/doc/UFsparse/UFconfig.txt
mv UFconfig/README.txt ../inst/doc/UFsparse/UFconfig.txt
  ## touch the file UFconfig/UFconfig.mk.  We use other configuration
  ## environment variables but this name is embedded in some Makefiles
touch UFconfig/UFconfig.mk
  ## install COLAMD/Source and COLAMD/Include directories

## 2) COLAMD -----------------------------------------------
tar zxf COLAMD.tar.gz COLAMD/Source/ COLAMD/Include/
Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("COLAMD/Source/Makefile")'
  ## install documentation for COLAMD
tar zxf COLAMD.tar.gz COLAMD/README.txt COLAMD/ChangeLog
mv COLAMD/README.txt ../inst/doc/UFsparse/COLAMD.txt
  ## install AMD/Source and AMD/Include directories

## 3) AMD --------------------------------------------------
tar zxf AMD.tar.gz AMD/Source AMD/Include AMD/README.txt
  ## restore the AMD/Source/Makefile
svn revert AMD/Source/Makefile
  ## install AMD documentation
mv AMD/README.txt ../inst/doc/UFsparse/AMD.txt
  ## remove Fortran source files and GNUMakefile
rm AMD/Source/*.f AMD/Source/GNUmakefile

## 4) CHOLMOD ----------------------------------------------
  ## install CHOLMOD source files
for d in Check Cholesky Core Include Lib MatrixOps Modify Partition Supernodal
do
    tar zxf ./CHOLMOD.tar.gz CHOLMOD/$d
done
  ## install CHOLMOD documentation in ../inst/doc/UFsparse
tar zxf ./CHOLMOD.tar.gz CHOLMOD/README.txt
mv CHOLMOD/README.txt ../inst/doc/UFsparse/CHOLMOD.txt

cp -p CHOLMOD/Lib/Makefile CHOLMOD/Lib/Makefile_CHOLMOD
Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("CHOLMOD/Lib/Makefile")'
mv    CHOLMOD/Lib/Makefile CHOLMOD/Lib/Makefile_pre
svn revert CHOLMOD/Lib/Makefile
  ##
ls -l CHOLMOD/Lib/Makefile_pre
echo 'now diff CHOLMOD/Lib/Makefile with CHOLMOD/Lib/Makefile_pre'
echo ' make changes as necessary, and then (later)'
echo ' rm CHOLMOD/Lib/Makefile_*' ; echo

## 5) SPQR -------------------------------------------------
  ## install SPQR source files
for d in Source Include Lib
do
    tar zxf ./SPQR.tar.gz SPQR/$d
done
  ## install CHOLMOD documentation in ../inst/doc/UFsparse
tar zxf ./SPQR.tar.gz SPQR/README.txt
mv SPQR/README.txt ../inst/doc/UFsparse/SPQR.txt

cp -p SPQR/Lib/Makefile SPQR/Lib/Makefile_SPQR
Rscript --vanilla -e 'source("scripts/fixup-fn.R")' -e 'fixup("SPQR/Lib/Makefile")'
mv    SPQR/Lib/Makefile SPQR/Lib/Makefile_pre
svn revert SPQR/Lib/Makefile
  ##
ls -l SPQR/Lib/Makefile_pre
echo 'now diff SPQR/Lib/Makefile with SPQR/Lib/Makefile_pre'
echo ' make changes as necessary, and then (later)'
echo ' rm SPQR/Lib/Makefile_*' ; echo

## ----- remove the downloaded tar files -------------------
rm CHOLMOD.tar.gz AMD.tar.gz COLAMD.tar.gz UFconfig.tar.gz SPQR.tar.gz
