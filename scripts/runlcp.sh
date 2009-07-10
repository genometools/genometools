#!/bin/sh

# set -e -x

if test $# -eq 0
then
  filenames="`find testdata -name '*.fna'` testdata/at1MB \
             `find testdata -name '*.fsa'`"
else
  if test $1 == 'valgrind'
  then
    RUNNER=valgrind.sh
    shift
  else
    RUNNER=time
  fi
  if test $# -eq 0
  then
    filenames="`find testdata -name '*.fna'` testdata/at1MB \
               `find testdata -name '*.fsa'`"
  else
    filenames=$*
  fi
fi

suffixeratornoidxnolcp()
{
  ${RUNNER} gt suffixerator -showtime -tis -suf -des -ssp -db ${filename} $*
}

suffixerator()
{
  suffixeratornoidxnolcp -lcp -indexname sfx-idx $*
}

suffixeratornolcp()
{
  suffixeratornoidxnolcp -indexname sfx-idx $*
}

suffixeratornoidx()
{
  suffixeratornoidxnolcp -lcp $*
}

sfxmap()
{
  gt dev sfxmap -lcp -suf $*
}

sfxmaponlysuf()
{
  gt dev sfxmap -suf $*
}

smallfiles="testdata/Random-Small.fna testdata/TTT-small.fna"

checksmallfile()
{
  filename=$1
  ret=0
  for cfc in ${smallfiles}
  do
    if test ${cfc} == ${filename}
    then
      ret=1
    fi
  done
  echo ${ret}
}

for filename in ${filenames}
do
  suffixerator ""
  sfxmap sfx-idx
  suffixerator -dir rev
  sfxmap sfx-idx
  suffixerator -maxdepth
  sfxmap sfx-idx
  maxdepth=`grep '^prefixlength=' sfx-idx.prj | sed -e 's/prefixlength=//'`
  maxdepth=`expr ${maxdepth} \* 2`
  suffixeratornoidx -maxdepth ${maxdepth} -indexname sfx-idx${maxdepth}
  sfxmap sfx-idx${maxdepth}
  suffixerator -parts 3
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth he
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth abs
  sfxmap sfx-idx
  gt dev guessprot ${filename}
  if test $? -eq 0
  then
    strandlist="fwd rev cpl rcl"
  else
    strandlist="fwd rev"
  fi
  for dir in ${strandlist}
  do
    for dc in 8 16 32
    do
      suffixeratornolcp -dir ${dir} -cmpcharbychar -dc ${dc}
      sfxmaponlysuf sfx-idx
      suffixeratornolcp -dir ${dir} -dc ${dc}
      sfxmaponlysuf sfx-idx
    done
  done
  xx=`checksmallfile ${filename}`
  if test ${xx} == '0'
  then
    suffixeratornolcp -pl 2
    sfxmaponlysuf sfx-idx
  fi
  # ${RUNNER} gt suffixerator -v -showtime -smap Transab -tis -suf -dc 64 -db testdata/fib25.fas.gz -indexname sfx-idx
  rm -f sfx-idx.* sfx-idx${maxdepth}.* 
done
echo "${filenames}"
