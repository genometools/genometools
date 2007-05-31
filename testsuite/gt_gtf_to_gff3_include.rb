Name "gt gtf2gff3 test"
Keywords "gt_gtf2gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test.gtf"
  run "diff #{$last_stdout} #{$testdata}gt_gtf_to_gff3_test.gff3"
end
