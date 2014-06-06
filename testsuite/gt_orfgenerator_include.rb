Name "gt orfgenerator missing input sequence"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator", :retval => 1, :maxtime => 120
  grep(last_stderr, /missing argument/)
end

Name "gt orfgenerator missing value for -min"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator -min foo foo.gff3", :retval => 1
  grep(last_stderr, /error/)
end

Name "gt orfgenerator missing value for -max"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator -max foo foo.gff3", :retval => 1
  grep(last_stderr, /error/)
end

Name "gt orfgenerator (longest for stop codon)"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator #{$testdata}/Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.orfs30.gff3"
end

Name "gt orfgenerator (all ORFs)"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator -all #{$testdata}/Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.orfs30.all.gff3"
end

Name "gt orfgenerator (empty sequence)"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator #{$testdata}/empty_seq.fas", :retval => 1
end

Name "gt orfgenerator (empty sequence in between)"
Keywords "gt_orfgenerator"
Test do
  run_test "#{$bin}gt dev orfgenerator #{$testdata}/second_empty.fas"
end