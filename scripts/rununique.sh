#!/usr/bin/env bash

if test $# -lt 4
then
  echo "Usage: $0 <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

function cerr() 
{
  $1
  if test $? -ne 0
  then
    echo "failure: ${1}"
    exit 1
  fi
}

function mkfmindex() 
{
  indexname=$1
  shift
  iiargs=$*
  cerr "../bin/gt mkfmindex -noindexpos -fmout ${indexname} -ii ${iiargs}"
}

function plain() 
{
  cerr "../bin/gt suffixerator -plain -tis -pl 1 -indexname $1 -smap $2 -db $3"
}

function suffixerator
{
  cerr "../bin/gt suffixerator $*"
}

function uniquesub
{
  cerr "../bin/gt uniquesub -output sequence querypos $*"
}

IDIR=../indexdir
mkdir -p ${IDIR}

minval=$1
maxval=$2
queryfile=$3
shift
shift
shift
rfiles=$*

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
  if test ! -f ${IDIR}/mkv-single.prj ||
     test ${rfiles} -nt ${IDIR}/mkv-single.prj   
  then
    suffixerator -indexname ${IDIR}/mkv-single -db ${rfiles} -bwt -pl 8
  else
    echo "# ${IDIR}/mkv-single is up to date"
  fi
  if test ! -f ${IDIR}/fmidx.fma ||
     test ${IDIR}/mkv-single.prj -nt ${IDIR}/fmidx.fma
  then
    mkfmindex ${IDIR}/fmidx ${IDIR}/mkv-single
    plain ${IDIR}/fmidx ${IDIR}/mkv-single.al1 ${IDIR}/mkv-single.bwt
  else
    echo "# ${IDIR}/fmidx is up to date"
  fi
else
  runmkfm=0
  fn=0
  for rfile in ${rfiles}
  do
    if test ! -f ${IDIR}/midx${fn}.prj || 
       test ${rfile} -nt ${IDIR}/midx${fn}.prj
    then
      suffixerator -indexname ${IDIR}/midx${fn} -db ${rfile} -suf -lcp -tis -pl 1
      runmkfm=1
    else
      echo "# ${IDIR}/midx${fn} is up to date"
    fi
    indexlist="${indexlist} ${IDIR}/midx${fn}"
    fn=`expr ${fn} + 1`
  done
  if test ${runmkfm} -eq 1
  then
    mkfmindex ${IDIR}/fmidx ${indexlist}
    plain ${IDIR}/fmidx ${IDIR}/fmidx.al1 ${IDIR}/fmidx.bwt
  else
    echo "# ${IDIR}/fmidx is up to date"
  fi
fi
fmindexname=${IDIR}/fmidx
uniquesub -fmi ${fmindexname} -query ${queryfile} -min ${minval} -max ${maxval}
