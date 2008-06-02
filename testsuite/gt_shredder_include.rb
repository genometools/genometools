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
