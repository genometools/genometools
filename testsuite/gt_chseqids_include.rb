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
