#!/bin/sh
## SVN_MATRIX_DIR is the directory produced by something like
## svn co svn+ssh://mmaechler@svn.r-forge.r-project.org/svnroot/matrix Matrix
SVN_MATRIX_DIR=~/R/D/R-forge/Matrix
##
checkDir () {
    dir="$1"
    Nam="$2"
    if [ ! -d $dir ]
    then
	if [ x$Nam != x ] ; then echo -n "$Nam="; fi
	echo "'$dir' is not a valid directory"
    fi
}

checkDir $SVN_MATRIX_DIR SVN_MATRIX_DIR
d=$SVN_MATRIX_DIR/pkg/Matrix/inst ; checkDir $d ; cd $d

doxygen -u Doxyfile
doxygen
d=$SVN_MATRIX_DIR/www/doxygen; checkDir $d ; cd $d
find . -path ./.svn -prune -o \( -mtime +1 -a -type f -exec rm {} \; \)
svn add *

cd $SVN_MATRIX_DIR
echo -n "svn cleanup in `pwd` : "
svn cleanup .
echo '[Ok]
Now trying to commit all : '
set -v
svn ci -m'after "doxygen -u" and "doxygen" and cleanup of www directory' pkg/Matrix/inst/Doxyfile www/doxygen
