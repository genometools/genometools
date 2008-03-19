Name "gt sequniq test "
Keywords "gt_sequniq"
Test do
  run_test "#{$bin}gt sequniq #{$testdata}foofoo.fas"
  run "diff #{$last_stdout} #{$testdata}foo.fas"
end
