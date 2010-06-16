Name "seqiterator (FASTA)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.fna"
  grep($last_stdout, /21 sequences/)
end

Name "seqiterator (GenBank)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.gbk"
  grep($last_stdout, /21 sequences/)
end

Name "seqiterator (EMBL)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}Atinsert.embl"
  grep($last_stdout, /21 sequences/)
end

Name "seqiterator distlen"
Keywords "gt_seqiterator"
Test do
  {10 => "860--869 1", 
   100 => "800--899 4",
   1000 => "0--999 21"}.each_pair do |bucketsize, expected|
    run_test "#{$bin}gt dev seqiterator -v -distlen"+
                      " -bucketsize #{bucketsize} #{$testdata}Atinsert.gbk"
    grep($last_stdout, expected)
  end
end

Name "seqiterator contigs"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -contigs #{$testdata}at1MB"
  grep($last_stdout, "number.*1952")
  grep($last_stdout, "total.*770425")
  grep($last_stdout, "average.*394")
  grep($last_stdout, "longest.*1102")
  grep($last_stdout, "N50.*421")
  grep($last_stdout, "smallest.*56")
end

Name "seqiterator fail (unknown file type)"
Keywords "gt_seqiterator"
Test do
  run_test "#{$bin}gt dev seqiterator -v -distlen #{$testdata}empty_file", \
           :retval => 1
  grep($last_stderr, /cannot guess file type/)
end
