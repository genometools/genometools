#!/bin/sh

cerr()
{
  $*
  if [ $? -ne 0 ]
  then
    echo "failure: $*"
    exit 1
  fi
}

for inputfile in `ls testdata/*.fna` ${AT} ${U8} ${ATK} `ls ${GTTESTDATA}/DNA-mix/Grumbach.fna/*.fna`
do
  if test ${inputfile} = 'testdata/TTT-small.fna'
  then
    echo "skip ${inputfile}"
  else
    if test ${inputfile} = 'testdata/Random-Small.fna'
    then
      echo "skip ${inputfile}"
    else
      echo "${inputfile}"
      cerr "scripts/checktallymer.sh ${inputfile}"
    fi
  fi
done

cd testsuite

env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb -keywords 'gt_greedyfwdmat'

cd ..
