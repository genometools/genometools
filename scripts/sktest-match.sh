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
  cerr "scripts/checktallymer.sh ${inputfile}"
done

cd testsuite

env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb -keywords 'gt_greedyfwdmat'

cd ..
