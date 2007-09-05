#!/bin/sh

if test $# -lt 4
then
  echo "Usage: $0 <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

fmdir=${PROJECT}/vstree/src/vstree/src/fmmatch

cerr() 
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

mkfmindex() 
{
  indexname=$1
  shift
  iiargs=$*
  cerr "${fmdir}/mkfmindex.x -noindexpos -fmout ${indexname} ${iiargs}"
}

mkvtree
{
  cerr "mkvtree.x $*"
}

uniquesub
{
  cerr "${fmdir}/uniquesub.x -output sequence querypos $*"
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
  if test ! -f ${IDIR}/mkv-single.rev.prj ||
     test ${rfiles} -nt ${IDIR}/mkv-single.rev.prj   
  then
    mkvtree -rev -dna -indexname ${IDIR}/mkv-single -db ${rfiles} -suf -tis -bwt -pl
  else
    echo "# ${IDIR}/mkv-single.rev is up to date"
  fi
  if test ! -f ${IDIR}/mkv-fm.fma ||
     test ${IDIR}/mkv-single.rev.prj -nt ${IDIR}/mkv-fm.fma
  then
    mkfmindex ${IDIR}/mkv-fm ${IDIR}/mkv-single.rev
  else
    echo "# ${IDIR}/mkv-fm is up to date"
  fi
else
  runmkfm=0
  fn=0
  for rfile in ${rfiles}
  do
    if test ! -f ${IDIR}/midx${fn}.prj || 
       test ${rfile} -nt ${IDIR}/midx${fn}.prj
    then
      mkvtree -rev -indexname ${IDIR}/midx${fn} -db ${rfile} -suf -lcp -tis -pl
      runmkfm=1
    else
      echo "# ${IDIR}/midx${fn} is up to date"
    fi
    indexlist="${indexlist} ${IDIR}/midx${fn}"
    fn=`expr ${fn} + 1`
  done
  if test ${runmkfm} -eq 1
  then
    mkfmindex ${IDIR}/mkv-fm ${indexlist}
  else
    echo "# ${IDIR}/mkv-fm is up to date"
  fi
fi
fmindexname=${IDIR}/mkv-fm
uniquesub -min ${minval} -max ${maxval} -smallindex ${fmindexname} ${queryfile}
