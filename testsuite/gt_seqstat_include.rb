Name "gt seqstat (FASTA)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.fna"
  grep($last_stdout, /21 sequences/)
end

Name "gt seqstat (GenBank)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.gbk"
  grep($last_stdout, /21 sequences/)
end

Name "gt seqstat (EMBL)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.embl"
  grep($last_stdout, /21 sequences/)
end

Name "gt seqstat distlen"
Keywords "gt_seqstat"
Test do
  {10 => "860--869 1", 
   100 => "800--899 4",
   1000 => "0--999 21"}.each_pair do |bucketsize, expected|
    run_test "#{$bin}gt seqstat -v -distlen"+
                      " -b #{bucketsize} #{$testdata}Atinsert.gbk"
    grep($last_stdout, expected)
  end
end

Name "gt seqstat contigs"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -contigs #{$testdata}at1MB"
  grep($last_stdout, "number.*1952")
  grep($last_stdout, "total.*770425")
  grep($last_stdout, "average.*394")
  grep($last_stdout, "longest.*1102")
  grep($last_stdout, "N50.*421")
  grep($last_stdout, "smallest.*56")
end

Name "gt seqstat fail (unknown file type)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}empty_file", \
           :retval => 1
  grep($last_stderr, /cannot guess file type/)
end
