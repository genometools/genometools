Name "gt mmapandread test 1"
Keywords "gt_mmapandread"
Test do
  run_test("#{$bin}gt mmapandread", :retval => 1)
end

Name "gt mmapandread test 2"
Keywords "gt_mmapandread"
Test do
  run_test "#{$bin}gt mmapandread #{$testdata}gt_mmapandread_test_1.in"
end

Name "gt mmapandread test 3"
Keywords "gt_mmapandread"
Test do
  run_test "#{$bin}gt mmapandread #{$testdata}gt_mmapandread_test_2.in"
end

Name "gt mmapandread test 4"
Keywords "gt_mmapandread"
Test do
  run_test "#{$bin}gt mmapandread #{$testdata}gt_mmapandread_test_3.in"
end

Name "gt mmapandread test 5"
Keywords "gt_mmapandread"
Test do
  run "ln -s /dev/null devnulllink" 
  run_test "#{$bin}gt mmapandread devnulllink"
end

Name "gt mmapandread test 6"
Keywords "gt_mmapandread"
Test do
  run_test "#{$bin}gt mmapandread #{$testdata}gt_mmapandread_test_1.in #{$testdata}gt_mmapandread_test_2.in #{$testdata}gt_mmapandread_test_3.in"
end
