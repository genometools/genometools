Name "gt chseqids -help"
Keywords "gt_chseqids"
Test do
  run_test "#{$bin}gt chseqids -help"
end

Name "gt chseqids -noop"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids -noop", :retval => 1)
end

Name "gt chseqids empty_file"
Keywords "gt_chseqids"
Test do
  run_test("#{$bin}gt chseqids #{$testdata}empty_file", :retval => 1)
  grep $last_stderr, "not defined"
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
