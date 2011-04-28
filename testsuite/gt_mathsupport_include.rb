Name "mathsupport Lua rand_max"
Keywords "gt_scripts mathsupport"
Test do
  run_test "#{$bin}gt -seed 12345 #{$testdata}/gtscripts/mathsupport.lua"
end

Name "mathsupport Lua rand_max (wrong seed)"
Keywords "gt_scripts mathsupport"
Test do
  run_test "#{$bin}gt -seed 54321 #{$testdata}/gtscripts/mathsupport.lua", \
           :retval => 1
end
