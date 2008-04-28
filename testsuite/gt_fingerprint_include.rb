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
