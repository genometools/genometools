#!/bin/sh

regularfiles=at1MB
repetitivefiles=mfd

code2file()
{
  case $1 in
    at1MB)
      echo "testdata/at1MB";;
    mfd)
      echo "${HOME}/seqcmpprojects/MouthFootDisease/mfdallpart.fna.gz";;
    *)
      echo "$0: illegal filecode $1"
      exit 1;;
  esac
}

checkregular()
{
  filename=$1
  for cfc in ${regularfiles}
  do
    if test ${cfc} = ${filename}
    then
      return 0
    fi
  done
  return 1
}

suffixerator()
{
  fc=$1
  filename=`code2file $1`
  shift
  printf "# RUN $fc $*\n"
  ${RUNNER} gt suffixerator -showtime -indexname sfx-id -dna -tis -suf -db ${filename} $* | egrep '# TIME overall|# space peak'
}

for rfc in $regularfiles $repetitivefiles
do
  checkregular ${rfc}
  if test $? -eq 0
  then
    suffixerator ${rfc} -cmpcharbychar ""
    suffixerator ${rfc} ""
  fi
  for dc in 8 32 128
  do
    suffixerator ${rfc} -cmpcharbychar -dc ${dc}
    suffixerator ${rfc} -dc ${dc}
  done
  rm -f sfx-idx.*
done
