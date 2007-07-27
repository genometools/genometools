#!/usr/bin/env bash

if test $# -lt 5
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
  if test $? -ne 0
  then
    echo "failure: ${1}"
    exit 1
  else
    echo "# success $1"
  fi
}

function mkfmindex() 
{
  indexname=$1
  shift
  iiargs=$*
  cerr "${gtpath} mkfmindex -noindexpos -fmout ${indexname} -ii ${iiargs}"
}

function plain() 
{
  cerr "${gtpath} suffixerator -plain -tis -indexname $1 -smap $2 -db $3"
}

function suffixerator
{
  cerr "${gtpath} suffixerator $*"
}

function uniquesub
{
  cerr "${gtpath} uniquesub -output sequence querypos $*"
}

IDIR=../indexdir
mkdir -p ${IDIR}

numofrfiles=0
for rfile in ${rfiles}
do
  if test -f ${rfile}
  then
    numofrfiles=`expr ${numofrfiles} + 1`
  else
    echo "# $0: file \"${rfile}\" does not exist"
    exit 1
  fi
done

if test ${numofrfiles} -eq 1
then
  if test ! -f ${IDIR}/sfx-single.rev.prj ||
     test ${rfiles} -nt ${IDIR}/sfx-single.rev.prj   
  then
    suffixerator -dir rev -dna -indexname ${IDIR}/sfx-single.rev -db ${rfiles} -bwt -pl
  else
    echo "# ${IDIR}/sfx-single.rev is up to date"
  fi
  if test ! -f ${IDIR}/sfx-fm.fma ||
     test ${IDIR}/sfx-single.rev.prj -nt ${IDIR}/sfx-fm.fma
  then
    mkfmindex ${IDIR}/sfx-fm ${IDIR}/sfx-single.rev
    plain ${IDIR}/sfx-fm ${IDIR}/sfx-single.rev.al1 ${IDIR}/sfx-single.rev.bwt
  else
    echo "# ${IDIR}/sfx-fm is up to date"
  fi
else
  runmkfm=0
  fn=0
  for rfile in ${rfiles}
  do
    if test ! -f ${IDIR}/midx${fn}.prj || 
       test ${rfile} -nt ${IDIR}/midx${fn}.prj
    then
      suffixerator -dna -indexname ${IDIR}/midx${fn} -db ${rfile} -suf -lcp -tis -pl
      runmkfm=1
    else
      echo "# ${IDIR}/midx${fn} is up to date"
    fi
    indexlist="${indexlist} ${IDIR}/midx${fn}"
    fn=`expr ${fn} + 1`
  done
  if test ${runmkfm} -eq 1
  then
    mkfmindex ${IDIR}/sfx-fm ${indexlist}
    plain ${IDIR}/sfx-fm ${IDIR}/sfx-fm.al1 ${IDIR}/sfx-fm.bwt
  else
    echo "# ${IDIR}/sfx-fm is up to date"
  fi
fi
fmindexname=${IDIR}/sfx-fm
uniquesub -fmi ${fmindexname} -query ${queryfile} -min ${minval} -max ${maxval}
