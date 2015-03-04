Name "gt interfeat test"
Keywords "gt_interfeat"
Test do
  run_test "#{$bin}gt interfeat #{$testdata}addintrons.gff3"
  run "diff #{last_stdout} #{$testdata}addintrons.out"
end

Name "gt interfeat test (multi-feature)"
Keywords "gt_interfeat"
Test do
  run_test "#{$bin}gt interfeat -outside EST_match -inter match_gap #{$testdata}interfeat_pseudo.gff3"
  run "diff #{last_stdout} #{$testdata}interfeat_pseudo.out"
end
