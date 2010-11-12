Name "gt md5s_to_seqids (U89959_sas.gff3)"
Keywords "gt_md5s_to_seqids"
Test do
  run_test "#{$bin}gt md5s_to_seqids -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests.fas -- " +
           "#{$testdata}U89959_sas.gff3md5"
  run "diff #{$last_stdout} #{$testdata}U89959_sas.gff3"
end

Name "gt md5s_to_seqids (U89959_csas.gff3)"
Keywords "gt_md5s_to_seqids"
Test do
  run_test "#{$bin}gt md5s_to_seqids -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests.fas -- " +
           "#{$testdata}U89959_csas.gff3md5"
  run "diff #{$last_stdout} #{$testdata}U89959_csas.gff3"
end
