Name "gt loccheck test (fail)"
Keywords "gt_loccheck"
Test do
  run_test "#{$bin}gt loccheck < #{$testdata}gt_loccheck_containment_fail.gff3"
  grep last_stderr, 'mRNA child range 1123458-1124606'
end

Name "gt loccheck test (ok)"
Keywords "gt_loccheck"
Test do
  run_test "#{$bin}gt loccheck < #{$testdata}standard_gene_as_dag.gff3"
end