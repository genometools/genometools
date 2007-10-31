Name "gt regioncov test 1 (no options)"
Keywords "gt_regioncov"
Test do
  run "#{$bin}gt dev regioncov #{$testdata}encode_known_genes_Mar07.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_regioncov_test_1.out"
end

Name "gt regioncov test 2 (-maxfeaturedist)"
Keywords "gt_regioncov"
Test do
  run "#{$bin}gt dev regioncov -maxfeaturedist 220000 #{$testdata}encode_known_genes_Mar07.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_regioncov_test_2.out"
end
