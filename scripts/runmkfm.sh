#!/bin/sh

if [ $# -lt 3 ]
then
  echo "Usage: $0 <bothdirections> <idir> <fmindex> <file1> [file2 file3 ..]"
  exit 1
fi

bothdirections=$1
idir=$2
fmindex=${idir}/$3
shift
shift
shift
rfiles=$*

cerr() 
{
  $*
  if [ $? -ne 0 ]
  then
    echo "failure: $*"
    exit 1
  fi
}

suffixerator()
{
  cerr "../bin/gt suffixerator $*"
}

makesuftab()
{
  if [ $1 = 'rev' ]
  then
    suffixerator -dna -lcp -tis -suf -pl -dir rev -indexname ${idir}/$2.rev -db $3
  else
    suffixerator -dna -lcp -tis -suf -pl -dir cpl -indexname ${idir}/$2.cpl -db $3
  fi
}

plain() 
{
  cerr "../bin/gt suffixerator -plain -tis -indexname $1 -smap $1.al1 -db $1.bwt"
}

mkfmindex() 
{
  indexname=$1
  shift
  iiargs=$*
  cerr "../bin/gt mkfmindex -size small -noindexpos -fmout ${indexname} -ii ${iiargs}"
}

needsrebuild=0

needconstruction1()
{
  if [ ! -f ${idir}/$2 ] ||
     [ $1 -nt ${idir}/$2 ]
  then
    $3
    needsrebuild=1
  else
    echo "# ${idir}/$2 is up to date"
  fi
}

mkdir -p ${idir}

indexlist=""
for rfile in ${rfiles}
do
  indexname=`basename ${rfile}`
  needconstruction1 ${rfile} \
                    "${indexname}.rev.prj" \
                    "makesuftab rev ${indexname} ${rfile}"
  indexlist="${indexlist} ${idir}/${indexname}.rev"
  if [ $bothdirections -eq 1 ]
  then
    needconstruction1 ${rfile} \
                      "${indexname}.cpl.prj" \
                      "makesuftab cpl ${indexname} ${rfile}"
    indexlist="${indexlist} ${idir}/${indexname}.cpl"
  fi
done

if [ $needsrebuild -eq 1 ] ||
   [ ! -f ${fmindex}.fma ] ||
   [ ! -f ${fmindex}.fmd ] ||
   [ ! -f ${fmindex}.bwt ]
then
  mkfmindex ${fmindex} ${indexlist}
  needsrebuild=1
else
  echo "# ${fmindex}.fma is up to date"
  echo "# ${fmindex}.fmd is up to date"
  echo "# ${fmindex}.bwt is up to date"
fi

if [ $needsrebuild -eq 1 ]
then
  needconstruction1 ${fmindex}.bwt ${fmindex}.esq "plain ${fmindex}"
  needsrebuild=1
else
  echo "# ${fmindex}.esq is up to date"
fi
