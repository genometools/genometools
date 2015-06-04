#!/bin/sh

set -e -x

for filename in `cat filelist.txt`
do
  scripts/cmp-seedextend.rb ${filename} 20 20 10 30
done

/Users/stefan/genometools/testdata/Arabidopsis-C99826.fna
/Users/stefan/genometools/testdata/Atinsert.fna
/Users/stefan/genometools/testdata/Duplicate.fna
/Users/stefan/genometools/testdata/Ecoli-section1.fna
/Users/stefan/genometools/testdata/Ecoli-section2.fna
/Users/stefan/genometools/testdata/example_1.fa
/Users/stefan/genometools/testdata/nowildcardatend.fna
/Users/stefan/genometools/testdata/nowildcardatend_rev.fna
/Users/stefan/genometools/testdata/rcr_testseq.fa
/Users/stefan/genometools/testdata/Reads1.fna
/Users/stefan/genometools/testdata/Reads2.fna
/Users/stefan/genometools/testdata/Reads3.fna
/Users/stefan/genometools/testdata/Repfind-example.fna
/Users/stefan/genometools/testdata/sain.fna
/Users/stefan/genometools/testdata/Scaffold_102.fa
/Users/stefan/genometools/testdata/Small.fna
/Users/stefan/genometools/testdata/Smalldup.fna
/Users/stefan/genometools/testdata/test1.fasta
/Users/stefan/genometools/testdata/trna_glutamine.fna

/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/chntxx.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/hs5hcmvcg.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/humdystrop.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/humghcsa.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/humhbb.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/humhdabcd.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/humhprtb.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/mipacga.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/mpocpcg.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/mpomtcg.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/vaccg.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/Wildcards.fna
/Users/stefan/gttestdata/DNA-mix/Grumbach.fna/ychrIII.fna
