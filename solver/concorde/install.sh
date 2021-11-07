#!/usr/bin/env bash

# Download and compile Concorde

ARCHIVE=co031219.tgz
DIR=concorde
THIS_SCRIPT_DIR=`dirname "$0"`
CURRENT_DIR=`pwd`

cd $THIS_SCRIPT_DIR
if [ ! -f $ARCHIVE ] ; then
    wget http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/$ARCHIVE
fi

if [ ! -f $DIR ] ; then
    tar xvfz $ARCHIVE
fi

cd $DIR
cp ../lpcplex8.c LP
cp ../configure .
./configure --with-cplex=$1 --enable-pthreads 2>&1
make 2>&1
cd $CURRENT_DIR
