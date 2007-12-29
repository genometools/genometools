#!/bin/sh

# set -e -x

if [ $# -lt 4 ]
then
  echo "Usage: $0 <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

minval=$1
maxval=$2
queryfile=$3
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
  cerr "../bin/gt uniquesub -output sequence querypos $*"
}

idir=../indexdir
fmindex=Combined.fm

cerr "../scripts/runmkfm.sh ${idir} ${fmindex} ${rfiles}"

cerr "uniquesub -min ${minval} -max ${maxval} -fmi ${idir}/${fmindex} -query ${queryfile}"
