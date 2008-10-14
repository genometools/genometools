#!/bin/sh

checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
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
outoptions="-counts -pl -mersize ${mersize} -minocc 10 -maxocc 30"

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
${cmd} > tmp2
checkerror
cmd="valgrind.sh bin/gt tallymer mkindex ${outoptions} -indexname tyr-index sfxidx"
${cmd}
checkerror
cmd="valgrind.sh bin/gt tallymer search -strand fp -output qseqnum qpos counts sequence -test tyr-index ${AT}"
${cmd}
checkerror
rm -f tmp[12]
