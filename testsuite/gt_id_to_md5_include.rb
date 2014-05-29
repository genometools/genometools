Name "gt id_to_md5 (U89959_sas.gff3)"
Keywords "gt_id_to_md5"
Test do
  run_test "#{$bin}gt id_to_md5 -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests_unique.fas " +
           "-matchdesc #{$testdata}U89959_sas.gff3"
  run "diff #{last_stdout} #{$testdata}U89959_sas.gff3md5"
end

Name "gt id_to_md5 (U89959_csas.gff3)"
Keywords "gt_id_to_md5"
Test do
  run_test "#{$bin}gt id_to_md5 -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests_unique.fas " +
           "-matchdesc #{$testdata}U89959_csas.gff3"
  run "diff #{last_stdout} #{$testdata}U89959_csas.gff3md5"
end

Name "gt id_to_md5 (failure)"
Keywords "gt_id_to_md5"
Test do
  run_test "#{$bin}gt id_to_md5 -seqfiles " +
           "#{$testdata}U89959_genomic.fas #{$testdata}U89959_ests.fas " +
           "-matchdesc #{$testdata}U89959_csas.gff3", :retval => 1
  grep(last_stderr, "could match more than one sequence")
end
