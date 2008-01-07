#!/bin/sh

# set -e -x

if [ $# -lt 5 ]
then
  echo "Usage: $0 <gtbin> <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

gtbin=$1
minval=$2
maxval=$3
queryfile=$4
shift
shift
shift
shift
rfiles=$*

cerr() 
{
  $1
  if [ $? -ne 0 ]
  then
    echo "failure: ${1}"
    exit 1
  else
    echo "# success $1"
  fi
}

uniquesub()
{
  cerr "${gtbin} uniquesub -output sequence querypos $*"
}

idir=../indexdir
fmindex=Combined.fm

cerr "../scripts/runmkfm.sh ../bin/gt 1 ${idir} ${fmindex} ${rfiles}"

cerr "uniquesub -min ${minval} -max ${maxval} -fmi ${idir}/${fmindex} -query ${queryfile}"
