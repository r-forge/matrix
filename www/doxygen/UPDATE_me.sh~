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
# Remove all the doxygen files that are old, as
# they have not been *replaced* by a new version :
find . -path ./.svn -prune -o \( -mtime +1 -a -type f -exec rm {} \; \)
# but do "save myself":
svn revert UPDATE_me.sh
svn add * > svn-add.log 2>&1

cd $SVN_MATRIX_DIR
echo -n "svn cleanup in `pwd` (takes a while): "
svn cleanup .
echo '[Ok]
Now trying to commit all   (takes another while!): '
set -v
svn ci -m'after "doxygen -u" and "doxygen" and cleanup of www directory' pkg/Matrix/inst/Doxyfile www/doxygen
