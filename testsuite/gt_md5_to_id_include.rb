require 'fileutils'

Name "gt md5_to_id (U89959_sas.gff3, old)"
Keywords "gt_md5_to_id"
Test do
  FileUtils.copy "#{$testdata}U89959_ests.fas", "."
  FileUtils.copy "#{$testdata}U89959_genomic.fas", "."
  run_test "#{$bin}gt md5_to_id -seqfiles " +
           "U89959_genomic.fas U89959_ests.fas -- " +
           "#{$testdata}U89959_sas.gff3md5old"
  run "diff #{last_stdout} #{$testdata}U89959_sas.gff3"
end

Name "gt md5_to_id (U89959_sas.gff3)"
Keywords "gt_md5_to_id"
Test do
  run_test "#{$bin}gt md5_to_id #{$testdata}U89959_sas.gff3md5"
  run "diff #{last_stdout} #{$testdata}U89959_sas.gff3"
end

Name "gt md5_to_id (U89959_csas.gff3, old)"
Keywords "gt_md5_to_id"
Test do
  FileUtils.copy "#{$testdata}U89959_ests.fas", "."
  FileUtils.copy "#{$testdata}U89959_genomic.fas", "."
  run_test "#{$bin}gt md5_to_id -seqfiles " +
           "U89959_genomic.fas U89959_ests.fas -- " +
           "#{$testdata}U89959_csas.gff3md5old"
  run "diff #{last_stdout} #{$testdata}U89959_csas.gff3"
end

Name "gt md5_to_id (U89959_csas.gff3)"
Keywords "gt_md5_to_id"
Test do
  run_test "#{$bin}gt md5_to_id #{$testdata}U89959_csas.gff3md5"
  run "diff #{last_stdout} #{$testdata}U89959_csas.gff3"
end
