Name "gt seqids_to_md5s (U89959_sas.gff3)"
Keywords "gt_seqids_to_md5s"
Test do
  run_test "#{$bin}gt dev seqids_to_md5s -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests.fas " +
           "-matchdesc #{$testdata}U89959_sas.gff3"
  run "diff #{$last_stdout} #{$testdata}U89959_sas.gff3md5"
end

Name "gt seqids_to_md5s (U89959_csas.gff3)"
Keywords "gt_seqids_to_md5s"
Test do
  run_test "#{$bin}gt dev seqids_to_md5s -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests.fas " +
           "-matchdesc #{$testdata}U89959_csas.gff3"
  run "diff #{$last_stdout} #{$testdata}U89959_csas.gff3md5"
end
