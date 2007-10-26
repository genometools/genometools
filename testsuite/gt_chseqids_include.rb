Name "gt chseqids -help"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt chseqids -noop"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

Name "gt chseqids empty mapping file"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids #{$testdata}empty_file", :retval => 1)
  grep $last_stderr, "not defined"
end

Name "gt chseqids empty gff3 file"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_1.chseqids #{$testdata}empty_file"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt chseqids corrupt.gff3"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids #{$testdata}gt_chseqids_test_1.chseqids #{$testdata}corrupt.gff3", :retval => 1)
end

Name "gt chseqids test 1"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_1.chseqids #{$testdata}gt_chseqids_test_1.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_1.out"
end

Name "gt chseqids test 2"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_2.chseqids #{$testdata}gt_chseqids_test_2.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_2.out"
end

Name "gt chseqids test 3"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids #{$testdata}gt_chseqids_test_3.chseqids #{$testdata}gt_chseqids_test_3.gff3", :retval => 1)
end

Name "gt chseqids test 4"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_4.chseqids #{$testdata}gt_chseqids_test_4.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_4.out"
end

Name "gt chseqids test 5"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_5.chseqids #{$testdata}gt_chseqids_test_5.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_5.out"
end

Name "gt chseqids test 5 (-sort)"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids -sort #{$testdata}gt_chseqids_test_5.chseqids #{$testdata}gt_chseqids_test_5.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_5.sorted_out"
end

Name "gt chseqids test 6"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids #{$testdata}gt_chseqids_test_6.chseqids #{$testdata}gt_chseqids_test_6.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_chseqids_test_6.out"
end
