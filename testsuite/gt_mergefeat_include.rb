Name "gt mergefeat test"
Keywords "gt_mergefeat mergefeat"
Test do
  run_test "#{$bin}gt mergefeat #{$testdata}mergefeat.gff3"
  run "diff #{last_stdout} #{$testdata}mergefeat.out"
end

Name "gt mergefeat test (no merge)"
Keywords "gt_mergefeat mergefeat"
Test do
  run_test "#{$bin}gt mergefeat #{$testdata}mergefeat_no_merge.gff3"
  run "diff #{last_stdout} #{$testdata}mergefeat_no_merge.gff3"
end
