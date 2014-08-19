Name "fingerprint"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint #{$testdata}U89959_ests.fas | sort | uniq"
  run "diff #{last_stdout} #{$testdata}U89959_ests.checklist_uniq"
end

Name "fingerprint (nonexistent file)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt fingerprint #{$testdata}nonexistent_file", :retval=> 1)
end

Name "fingerprint (case insensitive)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt fingerprint #{$testdata}U89959_ests_gi_8690080_soft_masked.fas")
  run "diff #{last_stdout} #{$testdata}U89959_ests_gi_8690080_unmasked.checklist"
end

Name "fingerprint -check (success)"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -check #{$testdata}U89959_ests.checklist " +
           "#{$testdata}U89959_ests.fas"
end

Name "fingerprint -check (stdin)"
Keywords "gt_fingerprint"
Test do
  run "cat #{$testdata}U89959_ests.checklist | #{$memcheck} #{$bin}gt " +
      "fingerprint -check - #{$testdata}U89959_ests.fas"
end

Name "fingerprint -check (failure)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt fingerprint -check " +
           "#{$testdata}U89959_ests.checklist_uniq " +
           "#{$testdata}U89959_ests.fas", :retval => 1)
  grep last_stderr, /fingerprint comparison failed/
end

Name "fingerprint -check (failure)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt sequniq #{$testdata}U89959_ests.fas | #{$memcheck}" +
           "#{$bin}gt fingerprint -check #{$testdata}U89959_ests.checklist -",
           :retval => 1)
  grep last_stderr, /fingerprint comparison failed/
end

Name "fingerprint -collisions"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -collisions #{$testdata}U89959_ests.fas"
end

Name "fingerprint -duplicates (found)"
Keywords "gt_fingerprint"
Test do
  run_test("#{$bin}gt fingerprint -duplicates #{$testdata}U89959_ests.fas",
           :retval => 1)
  grep last_stderr, /duplicates found/
end

Name "fingerprint -duplicates (none found)"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -duplicates #{$testdata}U89959_genomic.fas"
end

Name "fingerprint -extract (found)"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -extract 6d3b4b9db4531cda588528f2c69c0a57 " +
           "#{$testdata}U89959_ests.fas"
  run "diff #{last_stdout} #{$testdata}gt_fingerprint_extract.out"
end

Name "fingerprint -extract (not found)"
Keywords "gt_fingerprint"
Test do
  run_test "#{$bin}gt fingerprint -extract 7d3b4b9db4531cda588528f2c69c0a57 " +
           "#{$testdata}U89959_ests.fas", :retval => 1
  grep last_stderr, /could not find sequence with fingerprint/
end
