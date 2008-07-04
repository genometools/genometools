Name "gt mutate test 1"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt -seed 1 mutate #{$testdata}gt_mutate_test_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_mutate_test_1.out"
end

Name "gt mutate test 2"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt -seed 987654321 mutate #{$testdata}gt_mutate_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_mutate_test_2.out"
end

Name "gt mutate test both"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt -seed 555555555 mutate -rate 2 #{$testdata}gt_mutate_test_1.fas #{$testdata}gt_mutate_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_mutate_test_both.out"
end
