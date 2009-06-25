#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 [small|all]"
  exit 1
fi

repetitivefiles="mfd paradoxus"

case $1 in
  small) allfiles="at1MB ecoli1 ecoli2 yeast mfd"
         ;;
  all)   allfiles="at1MB ecoli1 ecoli2 yeast mfd dmel human2 human22 paradoxus"
         ;;
  *)     allfiles=$1
         ;;
esac

code2file()
{
  case $1 in
    at1MB)
      echo "testdata/at1MB";;
    mfd)
      echo "${HOME}/seqcmpprojects/MouthFootDisease/mfdallpart.fna.gz";;
    yeast)
      echo "${HOME}/seqcmpprojects/yeast/yeast.fna.gz";;
    human2)
      echo "${PROJECT}/biodata/genomes/H_sapiens-build36.54-2009/Homo_sapiens.NCBI36.54.dna.chromosome.02.fa.gz";;
    human22)
      echo "${PROJECT}/biodata/genomes/H_sapiens-build36.54-2009/Homo_sapiens.NCBI36.54.dna.chromosome.22.fa.gz";;
    dmel)
      echo "${HOME}/seqcmpprojects/d_mel/d_mel.fna.gz";;
    ecoli1)
      echo "${PROJECT}/biodata/genomes/Bacteria/Ecoli/ecoli.fna";;
    ecoli2)
      echo "${PROJECT}/biodata/genomes/Bacteria/Ecoli_O157_H7/AE005174.fna";;
    swiss1MB)
      echo "${GTTESTDATA}/swissprot/swiss1MB";;
    paradoxus)
      echo "/local/kurtz/sfx-test/data/S-paradoxus.fna";;
    *)
      echo "$0: illegal filecode $1"
      exit 1;;
  esac
}

checkrepetitive()
{
  filename=$1
  for cfc in ${repetitivefiles}
  do
    if test ${cfc} = ${filename}
    then
      return 1
    fi
  done
  return 0
}

suffixerator()
{
  fc=$1
  filename=`code2file $1`
  shift
  printf "# RUN $fc $*\n"
  ${RUNNER} gt suffixerator -v -showtime -indexname sfx-id -tis -suf -db ${filename} $* | egrep '# TIME overall|# space peak'
}

mkesa()
{
  fc=$1
  printf "# RUN $fc mkesa\n"
  filename=`code2file $1`
  runmkesa-sfx.sh ${filename}
}

for rfc in $allfiles
do
  fn=`code2file ${rfc}`
  if test ! -f ${fn}
  then
    echo "FAILURE: ${fn} does not exist"
    exit 1
  fi
done

# suffixerator ecoli2 -sat uint32 -dc 128
# exit 0

echo "# DATE `date +%Y-%m-%d-%H:%M`"
for rfc in $allfiles
do
  checkrepetitive ${rfc}
  if test $? -eq 0
  then
    suffixerator ${rfc} -cmpcharbychar ""
    suffixerator ${rfc} ""
  fi
  for dc in 32 128 256
  do
    suffixerator ${rfc} -cmpcharbychar -dc ${dc}
    suffixerator ${rfc} -dc ${dc}
  done
  mkesa ${rfc}
  rm -f sfx-idx.* mkesa-idx.*
done
