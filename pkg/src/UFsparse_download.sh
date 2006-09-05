#!/bin/sh

wget http://www.cise.ufl.edu/research/sparse/camd/current/CAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/amd/current/AMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/cholmod/current/CHOLMOD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/colamd/current/COLAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/ccolamd/current/CCOLAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/UFconfig/current/UFconfig.tar.gz
  ## install UFconfig.h file (now needed by some UFsparse libraries)
tar zxf UFconfig.tar.gz UFconfig/UFconfig.h UFconfig/README.txt
  ## Move the UFconfig/README.txt file to ../inst/doc/UFsparse/UFconfig.txt
mv UFconfig/README.txt ../inst/doc/UFsparse/UFconfig.txt
  ## touch the file UFconfig/UFconfig.mk.  We use other configuration
  ## environment variables but this name is embedded in some Makefiles
touch UFconfig/UFconfig.mk
  ## install source files for CCOLAMD
tar zxf CCOLAMD.tar.gz CCOLAMD/ccolamd.c CCOLAMD/ccolamd.h CCOLAMD/ccolamd_global.c
  ## install documentation for CCOLAMD
tar zxf CCOLAMD.tar.gz CCOLAMD/README.txt CCOLAMD/ChangeLog
mv CCOLAMD/README.txt ../inst/doc/UFsparse/CCOLAMD.txt
  ## install source files for COLAMD
tar zxf COLAMD.tar.gz COLAMD/colamd.c COLAMD/colamd.h COLAMD/colamd_global.c
  ## install documentation for COLAMD
tar zxf COLAMD.tar.gz COLAMD/README.txt COLAMD/ChangeLog
mv COLAMD/README.txt ../inst/doc/UFsparse/COLAMD.txt
  ## install AMD/Source and AMD/Include directories
tar zxf AMD.tar.gz AMD/Source AMD/Include AMD/README.txt AMD/Doc/AMD_UserGuide.pdf
  ## restore the AMD/Source/Makefile
svn revert AMD/Source/Makefile
  ## install AMD documentation
mv AMD/Doc/AMD_UserGuide.pdf ../inst/doc/UFsparse
rmdir AMD/Doc
mv AMD/README.txt ../inst/doc/UFsparse/AMD.txt
  ## remove Fortran source files and GNUMakefile
rm AMD/Source/*.f AMD/Source/GNUmakefile
  ## install CAMD/Source and CAMD/Include directories
tar zxf CAMD.tar.gz CAMD/Source CAMD/Include CAMD/README.txt CAMD/Doc/CAMD_UserGuide.pdf
  ## install CAMD documentation
mv CAMD/Doc/CAMD_UserGuide.pdf ../inst/doc/UFsparse
rmdir CAMD/Doc
mv CAMD/README.txt ../inst/doc/UFsparse/CAMD.txt
  ## remove Fortran source files and GNUMakefile
rm AMD/Source/*.f AMD/Source/GNUmakefile
  ## install CHOLMOD source files
for d in Check Cholesky Core Include Lib MatrixOps Modify Partition Supernodal
do
    tar zxf ./CHOLMOD.tar.gz CHOLMOD/$d
done
  ## install CHOLMOD documentation in ../inst/doc/UFsparse
tar zxf ./CHOLMOD.tar.gz CHOLMOD/Doc/UserGuide.pdf CHOLMOD/README.txt
mv CHOLMOD/Doc/UserGuide.pdf ../inst/doc/UFsparse/CHOLMOD_Userguide.pdf
rmdir CHOLMOD/Doc/
mv CHOLMOD/README.txt ../inst/doc/UFsparse/CHOLMOD.txt

  ## now svn diff the CHOLMOD/Lib/Makefile and the downloaded file and
  ## make changes as necessary

  ## remove the downloaded tar files
rm CHOLMOD.tar.gz AMD.tar.gz COLAMD.tar.gz CCOLAMD.tar.gz UFconfig.tar.gz
