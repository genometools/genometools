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
