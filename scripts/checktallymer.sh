#!/bin/sh

checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "${cmd}"
fi
}

cerr()
{
  $*
  if [ $? -ne 0 ]
  then
    echo "failure: $*"
    exit 1
  fi
}

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

mersize=20
inputfile=$1
outoptions="-counts -pl -mersize ${mersize} -minocc 2 -maxocc 30"
#PRECMD="valgrind.sh"
PRECMD="env -i"

cerr "bin/gt suffixerator -db ${inputfile} -tis -suf -lcp -pl -dna -indexname sfxidx"
cmd="env -i bin/gt tallymer mkindex -test -mersize ${mersize} sfxidx"
${cmd} > tmp1
checkerror
cerr "mkvtree.x -db ${inputfile} -tis -suf -lcp -pl -dna -indexname mkvidx"
cmd="tallymer-mkindex -mersize ${mersize} mkvidx" 
${cmd} > tmp2
checkerror
grep -v '^#' tmp2 > tmp2.2
mv tmp2.2 tmp2
cerr "cmp -s tmp1 tmp2"
rm -f tmp[12]
cmd="tallymer-mkindex ${outoptions} -indexname mkv-tyr-index mkvidx" 
${cmd}
checkerror
cmd="${PRECMD} bin/gt tallymer mkindex ${outoptions} -indexname tyr-index sfxidx"
${cmd}
checkerror

if test -s tyr-index.mct
then
  cmd="${PRECMD} bin/gt tallymer search -strand fp -output qseqnum qpos counts sequence -test tyr-index ${AT}"
  ${cmd} > tmp1
  checkerror
  cmd="tallymer-search -strand fp -output qseqnum qpos counts sequence mkv-tyr-index ${AT}"
  ${cmd} > tmp2
  checkerror
  grep -v '^#' tmp2 > tmp2.2
  mv tmp2.2 tmp2
  cerr "cmp -s tmp1 tmp2"
fi

rm -f tmp[12]
