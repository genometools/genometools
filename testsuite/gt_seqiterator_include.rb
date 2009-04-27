Name "seqiterator (FASTA)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.fna"
end

Name "seqiterator (GenBank)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.gbk"
end

Name "seqiterator (EMBL)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.embl"
end

Name "seqiterator fail (unknown file type)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}empty_file", \
           :retval => 1
  grep($last_stderr, /cannot guess file type/)
end

Name "seqiterator w/ qualities FASTQ (success)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads -showseq #{$testdata}test1.fastq"
end

Name "seqiterator w/ qualities FASTQ (non-FASTQ file)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}eden.gff3", \
           :retval => 1
  grep($last_stderr, /expected/)
end

Name "seqiterator w/ qualities FASTQ (invalid block start)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test2_wrong_begin.fastq", \
           :retval => 1
  grep($last_stderr, /expected/)
end

Name "seqiterator w/ qualities FASTQ (different seqnames)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test3_different_seqnames.fastq", \
           :retval => 1
  grep($last_stderr, "sequence description 'HWI-EAS306_9_FC305MP_6_1_1331" + \
                     "_1843' is not equal to qualities description 'HWI-EAS3" +\
                     "06_9_FC305MP_6_1_1331' in line")
end

Name "seqiterator w/ qualities FASTQ (different seqlengths 1)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test4_different_seqlengths.fastq", \
           :retval => 1
  grep($last_stderr, "lengths of character sequence and qualities sequence " + \
                     "differ")
end

Name "seqiterator w/ qualities FASTQ (different seqlengths 2)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test9_uneven_length.fastq", \
           :retval => 1
  grep($last_stderr, "qualities string of sequence length 33 is not ended " + \
                     "by newline")
end

Name "seqiterator w/ qualities FASTQ (tricky)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test5_tricky.fastq"
end

Name "seqiterator w/ qualities FASTQ (empty sequence)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test7_empty_seq.fastq", \
           :retval => 1
  grep($last_stderr, /empty sequence/)
end

Name "seqiterator w/ qualities FASTQ (premature end)"
Keywords "gt_seqiterator_qual"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test6_premature_end.fastq", \
           :retval => 1
  grep($last_stderr, /premature end/)
end
