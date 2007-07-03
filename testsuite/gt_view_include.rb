Name "gt view short test (with gff3_file_1_short.txt)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt view out.png < #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
  run "rm out.png"
end
