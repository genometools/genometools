Name "gt csa test"
Keywords "gt_csa"
Test do
  run_test("#{$bin}gt csa #{$testdata}/gt_csa_test_1", :retval => 1)
end

1.upto(4) do |i|
  Name "gt csa prob #{i}"
  Keywords "gt_csa"
  Test do
    run_test "#{$bin}gt csa #{$testdata}/gt_csa_prob_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}/gt_csa_prob_#{i}.out"
  end
end
