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
