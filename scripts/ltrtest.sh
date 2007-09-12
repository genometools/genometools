#!/bin/sh
#
# Copyright (c) 2007 David Ellinghaus <dellinghaus@stud.zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

function checkerror() 
{
  $1
  if test $? -ne 0
  then
    echo "failure: ${1}"
    exit 1
  fi
}

# the make call normally used for development
#cd ../testsuite
#./testsuite.rb -keywords gt_ltr
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)

#for (( i=1; i<=9 ; i++ ))
#do	
#  filenames[$i]="${LTR}/gt_ltr/s_cer_tab/chr0${i}_tab/chr0${i}.fsa"
#done

#for (( i=10; i<=16 ; i++ ))
#do	
#  filenames[$i]="${LTR}/gt_ltr/s_cer_tab/chr${i}_tab/chr${i}.fsa"
#done
#filenames[17]="${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa"

#for (( i=1; i<=16 ; i++ ))
#do
#  checkerror "../bin/gt ltrharvest -index ${filenames[${i}]} -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v"
#done

#checkerror "../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chr02_tab/chr02.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput"

#checkerror "../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v"

#checkerror "../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 0 -overlaps all -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v"

#../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 0 -overlaps all -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -mintsd 5 -maxtsd 20 -v -longoutput

#../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 0 -overlaps all -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -motif tgca -motifmis 0 -v -longoutput

../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 0 -overlaps all -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput -mintsd 5 -maxtsd 20 

#extractregion.sh 29641 29645 ${LTR}/gt_ltr/s_cer_tab/chr02_tab/chr02.19970727.fsa
#extractregion.sh 35598 35602 ${LTR}/gt_ltr/s_cer_tab/chr02_tab/chr02.19970727.fsa

#extractregion.sh 29572 29692 ${LTR}/gt_ltr/s_cer_tab/chr02_tab/chr02.19970727.fsa 
#extractregion.sh 35537 35657 ${LTR}/gt_ltr/s_cer_tab/chr02_tab/chr02.19970727.fsa 
