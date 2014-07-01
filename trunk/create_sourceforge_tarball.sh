#!/bin/bash
# This script produces the tarred, zipped files to upload to sourceforge site for the enumeration code
TEMPDIR=~/tarballtemp
rm -rf $TEMPDIR
mkdir $TEMPDIR
cp README $TEMPDIR
cd $TEMPDIR
svn co svn://hagrid.byu.edu/svn/enumlib/trunk enumlib/trunk
svn co svn://hagrid.byu.edu/svn/celib/trunk  celib/trunk
tar zcvf enum.tar.gz * 
