Name "gt seed_extend selfcompare"
Keywords "gt_seed_extend seed extend encseq"
Test do
  run_test "#{$bin}gt encseq encode -des no -sds no -md5 no #{$testdata}gt_bioseq_succ_3.fas"
  run_test "#{$bin}gt seed_extend -k 10 gt_bioseq_succ_3.fas"
  run "diff #{last_stdout} #{$testdata}seedextend1.out"
end

Name "gt seed_extend many reads k=14"
Keywords "gt_seed_extend encseq"
Test do
  run_test "#{$bin}gt encseq encode -des no -sds no -md5 no #{$testdata}condenser/varlen_longer_ids_200.fas"
  run_test "#{$bin}gt seed_extend -k 14 varlen_longer_ids_200.fas"
  run "diff #{last_stdout} #{$testdata}seedextend2.out"
end

Name "gt seed_extend mirror"
Keywords "gt_seed_extend encseq mirror"
Test do
  run_test "#{$bin}gt seed_extend -k 8 -mirror #{$testdata}foo.64"
  run "diff #{last_stdout} #{$testdata}seedextend3.out"
end

