Name "gt mutate test 1"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt mutate #{$testdata}gt_mutate_test_1.fas"
  grep($last_stdout, /mutated/);
end

Name "gt mutate test 2"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt mutate #{$testdata}gt_mutate_test_1.fas"
  grep($last_stdout, /mutated/);
end

Name "gt mutate test both"
Keywords "gt_mutate"
Test do
  run_test "#{$bin}gt mutate -rate 2 #{$testdata}gt_mutate_test_1.fas #{$testdata}gt_mutate_test_2.fas"
  grep($last_stdout, /mutated/);
end
