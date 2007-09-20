#!/bin/sh
#
# Copyright (c) 2007 David Ellinghaus <dellinghaus@stud.zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

# the make call normally used for development
#cd ../testsuite
#./testsuite.rb -keywords gt_ltr
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)
#cd ../

for (( i=1; i<=9 ; i++ ))
do	
  filenames[$i]="${LTR}/gt_ltr/s_cer_tab/chr0${i}_tab/chr0${i}.fsa"
done

for (( i=10; i<=16 ; i++ ))
do	
  filenames[$i]="${LTR}/gt_ltr/s_cer_tab/chr${i}_tab/chr${i}.fsa"
done

for (( i=1; i<=16 ; i++ ))
do
  ../bin/gt ltrharvest -index ${filenames[${i}]} -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v
done
../bin/gt ltrharvest -index ${LTR}/gt_ltr/s_cer_tab/chrAll_tab/chrAll_before-1997-10-01.fsa -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v

# test new output on yeast against old LTRharvest program output on yeast
diff /projects/gi/ltr/gt_ltr/results/s_cer/ellinghaus/Tab8_1/result-s_cer-All-Chr.fsa /projects/gi/ltr/results/s_cerevisiae/ellinghaus/Tab8_1/result-s_cer-All-Chr.fsa 

../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chr2L_tab/2L_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chr2R_tab/2R_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chr3L_tab/3L_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chr3R_tab/3R_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chr4_tab/4_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
../bin/gt ltrharvest -index /projects/gi/ltr/gt_ltr/d_mel_tab/chrX_tab/X_genomic_dmel_RELEASE3-1.FASTA -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -longoutput
