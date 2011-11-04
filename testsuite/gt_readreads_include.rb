Name "gt readreads (success)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads -showseq #{$testdata}test1.fastq"
end

Name "gt readreads colorspace (success)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads -colorspace -showseq "+
    "#{$testdata}solid_color_reads.fastq"
end

Name "gt readreads (non-FASTQ file)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}eden.gff3", \
           :retval => 1
  grep(last_stderr, /expected/)
end

Name "gt readreads (invalid block start)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test2_wrong_begin.fastq", \
           :retval => 1
  grep(last_stderr, /expected/)
end

Name "gt readreads (different seqnames)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test3_different_seqnames.fastq", \
           :retval => 1
  grep(last_stderr, "sequence description 'HWI-EAS306_9_FC305MP_6_1_1331" + \
                     "_1843' is not equal to qualities description 'HWI-EAS3" +\
                     "06_9_FC305MP_6_1_1331' in line")
end

Name "gt readreads (different seqlengths 1)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test4_different_seqlengths.fastq", \
           :retval => 1
  grep(last_stderr, "lengths of character sequence and qualities sequence " + \
                     "differ")
end

Name "gt readreads (different seqlengths 2)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test9_uneven_length.fastq", \
           :retval => 1
  grep(last_stderr, "qualities string of sequence length 33 is not ended " + \
                     "by newline")
end

Name "gt readreads (tricky)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads " + \
           "#{$testdata}test5_tricky.fastq"
end

Name "gt readreads (empty sequence)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test7_empty_seq.fastq", \
           :retval => 1
  grep(last_stderr, /empty sequence/)
end

Name "gt readreads (premature end)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test6_premature_end.fastq", \
           :retval => 1
  grep(last_stderr, /premature end/)
end

Name "gt readreads (multiline)"
Keywords "gt_readreads"
Test do
  run_test "#{$bin}gt dev readreads #{$testdata}test10_multiline.fastq"
end
