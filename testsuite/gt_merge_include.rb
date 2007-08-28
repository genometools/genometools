Name "gt merge test 1"
Keywords "gt_merge"
Test do
  run_test "#{$bin}gt merge #{$testdata}gff3_file_1_short.txt"
end

Name "gt merge test 2"
Keywords "gt_merge"
Test do
  run_test "#{$bin}gt merge #{$testdata}gff3_file_1_short.txt"
end

Name "gt merge test 3"
Keywords "gt_merge"
Test do
  run_test "#{$bin}gt merge #{$testdata}gt_merge_prob_1.in1 #{$testdata}gt_merge_prob_1.in2"
  run "diff #{$last_stdout} #{$testdata}gt_merge_prob_1.out"
end

Name "gt merge test 4"
Keywords "gt_merge"
Test do
  run_test "#{$bin}gt merge #{$testdata}gt_merge_prob_2.in1 #{$testdata}gt_merge_prob_2.in2"
  run "diff #{$last_stdout} #{$testdata}gt_merge_prob_2.out"
end

Name "gt merge unsorted file"
Keywords "gt_merge"
Test do
  run_test("#{$bin}gt merge #{$testdata}unsorted_gff3_file.txt", :retval => 1)
  grep($last_stderr, "is not sorted")
end
