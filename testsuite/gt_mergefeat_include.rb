Name "gt mergefeat test"
Keywords "gt_mergefeat mergefeat"
Test do
  run_test "#{$bin}gt mergefeat #{$testdata}mergefeat.gff3"
  run "diff #{$last_stdout} #{$testdata}mergefeat.out"
end
