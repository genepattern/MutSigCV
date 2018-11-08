#!/bin/sh
execDir=`dirname $0`

mkdir tmpClean1
mkdir tmpClean2
base1=`basename $1`
base2=`basename $2`
/usr/bin/perl $execDir/cleanStdout.pl $1 > tmpClean1/$base1
/usr/bin/perl $execDir/cleanStdout.pl $2 > tmpClean2/$base2

diff -q tmpClean1/$base1 tmpClean2/$base2
status=$?
rm -rf tmpClean1 tmpClean2

exit $status