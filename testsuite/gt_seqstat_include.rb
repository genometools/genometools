Name "gt seqstat (FASTA)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.fna"
  grep(last_stdout, /21 sequences/)
end

Name "gt seqstat (GenBank)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.gbk"
  grep(last_stdout, /21 sequences/)
end

Name "gt seqstat (EMBL)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}Atinsert.embl"
  grep(last_stdout, /21 sequences/)
end

Name "gt seqstat distlen"
Keywords "gt_seqstat"
Test do
  {10 => "860--869 1",
   100 => "800--899 4",
   1000 => "0--999 21"}.each_pair do |bucketsize, expected|
    run_test "#{$bin}gt seqstat -v -distlen"+
                      " -b #{bucketsize} #{$testdata}Atinsert.gbk"
    grep(last_stdout, expected)
  end
end

Name "gt seqstat -contigs"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -contigs #{$testdata}at1MB"
  grep(last_stdout, "number.*1952")
  grep(last_stdout, "total.*770425")
  grep(last_stdout, "mean.*394")
  grep(last_stdout, "median.*399")
  grep(last_stdout, "longest.*1102")
  grep(last_stdout, "shortest.*56")
  grep(last_stdout, "contigs > 500.*198")
  grep(last_stdout, "contigs > 1K.*3")
  grep(last_stdout, "contigs > 10K.*0")
  grep(last_stdout, "contigs > 100K.*0")
  grep(last_stdout, "contigs > 1M.*0")
  grep(last_stdout, "N50.*421")
  grep(last_stdout, "L50.*797")
end

Name "gt seqstat -contigs -genome"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -genome 770425 -contigs #{$testdata}at1MB"
  grep(last_stdout, "total.*770425")
  grep(last_stdout, "genome.*770425")
  grep(last_stdout, "genome.*100.00 %")
  grep(last_stdout, "N50.*421")
  grep(last_stdout, "NG50.*421")
  grep(last_stdout, "L50.*797")
  grep(last_stdout, "LG50.*797")
  run_test "#{$bin}gt seqstat -genome 100000 -contigs #{$testdata}at1MB"
  grep(last_stdout, "genome.*100000")
  grep(last_stdout, "genome.*770.42 %")
  grep(last_stdout, "N50.*421")
  grep(last_stdout, "NG50.*591")
  grep(last_stdout, "L50.*797")
  grep(last_stdout, "LG50.*74")
  run_test "#{$bin}gt seqstat -genome 1000000 -contigs #{$testdata}at1MB"
  grep(last_stdout, "genome.*1000000")
  grep(last_stdout, "genome.*77.04 %")
  grep(last_stdout, "N50.*421")
  grep(last_stdout, "NG50.*386")
  grep(last_stdout, "L50.*797")
  grep(last_stdout, "LG50.*1086")
  run_test "#{$bin}gt seqstat -genome 1600000 -contigs #{$testdata}at1MB"
  grep(last_stdout, "genome.*1600000")
  grep(last_stdout, "genome.*48.15 %")
  grep(last_stdout, "N50.*421")
  grep(last_stdout, "NG50.*n\.a\.")
  grep(last_stdout, "L50.*797")
  grep(last_stdout, "LG50.*n\.a\.")
end

Name "gt seqstat fail (unknown file type)"
Keywords "gt_seqstat"
Test do
  run_test "#{$bin}gt seqstat -v -distlen #{$testdata}empty_file", \
           :retval => 1
  grep(last_stderr, /cannot guess file type/)
end
