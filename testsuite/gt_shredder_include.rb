Name "gt shredder"
Keywords "gt_shredder"
Test do
  run_test "#{$bin}gt shredder #{$testdata}U89959_genomic.fas"
end

Name "gt shredder (-coverage)"
Keywords "gt_shredder"
Test do
  run_test "#{$bin}gt shredder -coverage 5 #{$testdata}U89959_genomic.fas"
end

Name "gt shredder (-overlap)"
Keywords "gt_shredder"
Test do
  run_test "#{$bin}gt shredder -overlap 500 #{$testdata}U89959_genomic.fas"
end

Name "gt shredder (nonexistent file)"
Keywords "gt_shredder"
Test do
  run_test("#{$bin}gt shredder #{$testdata}nonexistent_file", :retval => 1)
end

Name "gt shredder (wrong argument usage)"
Keywords "gt_shredder"
Test do
  run_test("#{$bin}gt shredder -maxlength 15  #{$testdata}U89959_genomic.fas",
           :retval => 1)
  grep $last_stderr, /-minlength must be <= than -maxlength/
end
