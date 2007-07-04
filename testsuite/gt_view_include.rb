Name "gt view short test"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt view short test (stdin)"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png < #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt view prob 1"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png #{$testdata}gt_view_prob_1.gff3"
  run "test -e out.png"
end
