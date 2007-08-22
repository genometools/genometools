Name "gt eval test 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval #{$testdata}/gt_eval_test_1.in #{$testdata}/gt_eval_test_1.in"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_test_1.out"
end

2.upto(8) do |i|
  Name "gt eval test #{i}"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval #{$testdata}/gt_eval_test_#{i}.reality #{$testdata}/gt_eval_test_#{i}.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_#{i}.out"
  end
end

Name "gt eval prob 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval #{$testdata}/gt_eval_prob_1.reality #{$testdata}/gt_eval_prob_1.prediction"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_prob_1.out"
end

Name "gt eval prob 1 (swapped)"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval #{$testdata}/gt_eval_prob_1.prediction #{$testdata}/gt_eval_prob_1.reality"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_prob_1.out_swapped"
end

Name "gt eval -ltr test 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval -ltr #{$testdata}/gt_eval_ltr_test_1.in #{$testdata}/gt_eval_ltr_test_1.in"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_ltr_test_1.out"
end

2.upto(9) do |i|
  Name "gt eval -ltr test #{i}"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval -ltr #{$testdata}/gt_eval_ltr_test_#{i}.reality #{$testdata}/gt_eval_ltr_test_#{i}.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_ltr_test_#{i}.out"
  end
end

