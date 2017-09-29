require "fileutils"
Name "gt seqtransform"
Keywords "gt_seqtransform (invariant)"
Test do
  FileUtils.copy "#{$testdata}nGASP/protein_100.fas", "."
  run_test "#{$bin}gt seqtransform protein_100.fas"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}nGASP/protein_100.fas"
end

Name "gt seqtransform"
Keywords "gt_seqtransform -addstopaminos"
Test do
  FileUtils.copy "#{$testdata}nGASP/protein_100.fas", "."
  run_test "#{$bin}gt seqtransform -addstopaminos protein_100.fas"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}nGASP/protein_100_with_stop.fas"
end

Name "gt seqtransform"
Keywords "gt_seqtransform -addstopaminos (invariant)"
Test do
  FileUtils.copy "#{$testdata}nGASP/protein_100_with_stop.fas", "."
  run_test "#{$bin}gt seqtransform -addstopaminos protein_100_with_stop.fas"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}nGASP/protein_100_with_stop.fas"
end
