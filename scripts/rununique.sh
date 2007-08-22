#!/usr/bin/env bash

# set -e -x

if [ $# -lt 5 ]
then
  echo "Usage: $0 <gtpath> <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

gtpath=$1
minval=$2
maxval=$3
queryfile=$4
shift
shift
shift
shift
rfiles=$*

function cerr() 
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

function uniquesub
{
  cerr "${gtpath} uniquesub -output sequence querypos $*"
}

idir=../indexdir
fmindex=Combined.fm

cerr "runmkfm.sh ${gtpath} ${idir} ${fmindex} ${rfiles}"

cerr "uniquesub -min ${minval} -max ${maxval} -fmi ${idir}/${fmindex} -query ${queryfile}"
