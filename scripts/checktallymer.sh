#!/bin/sh

checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "run ${cmd}"
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

cerr "bin/gt suffixerator -sat direct -db ${inputfile} -algbds 3 43 120 -tis -suf -lcp -pl -dna -indexname sfxidx"
cmd="env -i bin/gt tallymer mkindex -test -mersize ${mersize} -esa sfxidx"
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
cmd="${PRECMD} bin/gt tallymer mkindex ${outoptions} -indexname tyr-index -esa sfxidx"
${cmd}
checkerror

if test -s tyr-index.mct
then
  cmd="${PRECMD} bin/gt tallymer search -strand fp -output qseqnum qpos counts sequence -test -tyr tyr-index -q ${AT}"
  ${cmd} > tmp1
  checkerror
  cmd="tallymer-search -strand fp -output qseqnum qpos counts sequence mkv-tyr-index ${AT}"
  ${cmd} > tmp2
  checkerror
  grep -v '^#' tmp2 > tmp2.2
  mv tmp2.2 tmp2
  cerr "cmp -s tmp1 tmp2"
  if test ${inputfile} = "testdata/Duplicate.fna"
  then
    echo "skip ${inputfile}"
  else
    cmd="${PRECMD} bin/gt tallymer occratio -minmersize 10 -maxmersize 500 -esa sfxidx -output total unique nonunique nonuniquemulti relative"
    ${cmd} > tmp1
    checkerror
    grep -v '^#' tmp1 > tmp1.2
    mv tmp1.2 tmp1
    cmd="tallymer-occratio -minmersize 10 -maxmersize 500 -output total unique nonunique nonuniquemulti relative mkvidx"
    ${cmd} > tmp2
    checkerror
    grep -v '^#' tmp2 > tmp2.2
    mv tmp2.2 tmp2
    cerr "cmp -s tmp1 tmp2"
  fi
fi

rm -f tmp[12]
