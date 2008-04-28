Name "fingerprint -duplicates (found)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt fingerprint -duplicates #{$testdata}U89959_ests.fas",
           :retval => 1)
  grep $last_stderr, /duplicates found/
end

Name "fingerprint -duplicates (none found)"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -duplicates #{$testdata}U89959_genomic.fas"
end

Name "fingerprint -extract"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -extract 6d3b4b9db4531cda588528f2c69c0a57 " +
           "#{$testdata}U89959_ests.fas"
  run "diff #{$last_stdout} #{$testdata}gt_fingerprint_extract.out"
end
