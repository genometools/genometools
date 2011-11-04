Name "gt seqtransform"
Keywords "gt_seqtransform (invariant)"
Test do
  run_test "#{$bin}gt seqtransform #{$testdata}nGASP/protein_100.fas"
  run "diff #{last_stdout} #{$testdata}nGASP/protein_100.fas"
end

Name "gt seqtransform"
Keywords "gt_seqtransform -addstopaminos"
Test do
  run_test "#{$bin}gt seqtransform -addstopaminos " +
           "#{$testdata}nGASP/protein_100.fas"
  run "diff #{last_stdout} #{$testdata}nGASP/protein_100_with_stop.fas"
end

Name "gt seqtransform"
Keywords "gt_seqtransform -addstopaminos (invariant)"
Test do
  run_test "#{$bin}gt seqtransform -addstopaminos " +
           "#{$testdata}nGASP/protein_100_with_stop.fas"
  run "diff #{last_stdout} #{$testdata}nGASP/protein_100_with_stop.fas"
end
