Name "gt seed_extend small polysequence with mirror"
Keywords "gt_seed_extend seedpair kmer polysequence xdrop extend"
Test do
  run_test "#{$bin}gt encseq encode -des no -sds no -md5 no " +
           "-indexname #{$testdata}small_poly.fas #{$testdata}small_poly.fas"
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -debug-kmer " +
           "-debug-seedpair #{$testdata}small_poly.fas"
  run "cmp -s #{last_stdout} #{$testdata}seedextend1.out"
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -extendxdrop 97 -l 10 " +
           "-mincoverage 11 #{$testdata}small_poly.fas"
  run "cmp -s #{last_stdout} #{$testdata}seedextend3.out"
end

Name "gt seed_extend at1MB memlimit"
Keywords "gt_seed_extend seedpair at1MB memlimit maxfreq verbose"
Test do
  run_test "#{$bin}gt encseq encode -des no -sds no -md5 no " +
           "-indexname #{$testdata}at1MB #{$testdata}at1MB"
  run_test "#{$bin}gt seed_extend -verify -debug-seedpair -memlimit 10MB " +
           "#{$testdata}at1MB"
  run "cmp -s #{last_stdout} #{$testdata}seedextend2.out"
  run_test "#{$bin}gt seed_extend -v -maxfreq 5 #{$testdata}at1MB"
  run "grep -v '^# Found and sorted 582230 k-mers' #{last_stdout}"
  run "grep -v '^# Collected and sorted 68577 seed pairs' #{last_stdout}"
end

