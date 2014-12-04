require "fileutils"

Name "gt seqmutate test 1"
Keywords "gt_seqmutate"
Test do
  FileUtils.copy("#{$testdata}gt_mutate_test_1.fas", ".")
  run_test "#{$bin}gt -seed 1 seqmutate gt_mutate_test_1.fas"
  run "diff #{last_stdout} #{$testdata}gt_mutate_test_1.out"
end

Name "gt seqmutate test 2"
Keywords "gt_seqmutate"
Test do
  FileUtils.copy("#{$testdata}gt_mutate_test_2.fas", ".")
  run_test "#{$bin}gt -seed 987654321 seqmutate gt_mutate_test_2.fas"
  run "diff #{last_stdout} #{$testdata}gt_mutate_test_2.out"
end

Name "gt seqmutate test both"
Keywords "gt_seqmutate"
Test do
  FileUtils.copy("#{$testdata}gt_mutate_test_1.fas", ".")
  FileUtils.copy("#{$testdata}gt_mutate_test_2.fas", ".")
  run_test "#{$bin}gt -seed 555555555 mutate -rate 2 gt_mutate_test_1.fas "\
    "gt_mutate_test_2.fas"
  run "diff #{last_stdout} #{$testdata}gt_mutate_test_both.out"
end
