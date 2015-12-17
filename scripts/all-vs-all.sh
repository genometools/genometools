#!/bin/sh

set -e -x
# set GTDIR as path of genometools directory
excludelist=Verysmall.fna,TTTN.fna,description_test.fastq,Random160.fna,description_test2.fastq,Copysorttest.fna,TTT-small.fna,paired_reads_1N_2.p.fastq,Random.fna,wildcard1.fastq,wildcard2.fastq,fastq_problem2.fasta,fastq_problem.fasta,Reads3.fna,invalid_with_pos256_q31.fastq,Small.fna,sg_reads.fastq,Reads2.fna,hop_moderate.fastq,hop_aggressive.fastq

for reference in `${GTDIR}/scripts/findfasta.rb -n -e ${excludelist}`
do
  bin/gt encseq encode -indexname ref-index ${reference}
  for query in `${GTDIR}/scripts/findfasta.rb -n -e ${excludelist}`
  do
    bin/gt encseq encode -indexname query-index ${query}
    for mode in "" -seed-display
    do
      bin/gt seed_extend $mode -v -ii ref-index -maxfreq 20 -qii query-index > tmp.matches
      bin/gt dev show_seedext -a -f tmp.matches
      bin/gt dev show_seedext -e -f tmp.matches
      bin/gt dev show_seedext -a -f tmp.matches -sort
      bin/gt dev show_seedext -e -f tmp.matches -sort
      bin/gt dev show_seedext -f tmp.matches -sort | ascQorder.rb -c
    done
  done
done
