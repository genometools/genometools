Name "gt -help"
Keywords "gt"
Test do
  run_test "#{$bin}gt -help"
end

Name "gt -noop"
Keywords "gt"
Test do
  run_test("#{$bin}gt -noop", :retval => 1)
end
