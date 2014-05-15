Name "gt testcode_filter"
Keywords "gt_testcode_filter"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev testcode_filter -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.tcode.gff3"
end

Name "gt testcode_filter (all ORFs)"
Keywords "gt_testcode_filter"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev testcode_filter -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.tcode.all.gff3"
end

Name "gt testcode_filter (high threshold)"
Keywords "gt_testcode_filter"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev testcode_filter -threshold 33 -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.header_only.gff3"
end

Name "gt testcode_filter (empty sequence)"
Keywords "gt_testcode_filter"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run_test "#{$bin}gt dev testcode_filter -seqfile Lmaj1_start.fas -matchdesc #{$testdata}/empty_file", :retval => 1
end

Name "gt testcode_filter (invalid sequence)"
Keywords "gt_testcode_filter"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev testcode_filter -seqfile nonexist -matchdesc < #{last_stdout}", :retval => 1
end