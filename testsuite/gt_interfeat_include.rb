Name "gt interfeat test"
Keywords "gt_interfeat"
Test do
  run_test "#{$bin}gt interfeat #{$testdata}addintrons.gff3"
  run "diff #{last_stdout} #{$testdata}addintrons.out"
end
