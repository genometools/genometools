1.upto(14) do |i|
  Name "gt cds test #{i}"
  Keywords "gt_cds"
  Test do
    run_test "#{$bin}gt cds -minorflen 1 -startcodon yes " +
             "-seqfile #{$testdata}gt_cds_test_#{i}.fas " +
             "#{$testdata}gt_cds_test_#{i}.in"
    run "diff #{last_stdout} #{$testdata}gt_cds_test_#{i}.out"
  end
end

Name "gt cds error message"
Keywords "gt_cds"
Test do
  run "#{$bin}gt gff3 -offset 1000 #{$testdata}gt_cds_test_1.in | " +
      "#{$bin}gt cds -seqfile #{$testdata}gt_cds_test_1.fas -", :retval => 1
  grep last_stderr, "Has the sequence-region to sequence mapping been defined correctly"
end

1.upto(14) do |i|
  Name "gt cds test #{i} (-usedesc)"
  Keywords "gt_cds usedesc"
  Test do
    run_test "#{$bin}gt cds -minorflen 1 -startcodon yes -usedesc " +
             "-seqfile #{$testdata}gt_cds_test_#{i}.fas " +
             "#{$testdata}gt_cds_test_#{i}.in"
    run "diff #{last_stdout} #{$testdata}gt_cds_test_#{i}.out"
  end
end

Name "gt cds test (description range)"
Keywords "gt_cds usedesc"
Test do
  run_test "#{$bin}gt cds -minorflen 1 -usedesc -seqfile " +
           "#{$testdata}gt_cds_test_descrange.fas " +
           "#{$testdata}gt_cds_test_descrange.in"
  run "diff #{last_stdout} #{$testdata}gt_cds_test_descrange.out"
end

Name "gt cds test (multi description)"
Keywords "gt_cds usedesc"
Test do
  run_test "#{$bin}gt cds -minorflen 1 -usedesc -seqfile " +
           "#{$testdata}gt_cds_descrange_multi.fas " +
           "#{$testdata}gt_cds_descrange_multi.in"
  run "diff #{last_stdout} #{$testdata}gt_cds_descrange_multi.out"
end

Name "gt cds test (multi description fail 1)"
Keywords "gt_cds usedesc"
Test do
  run_test("#{$bin}gt cds -usedesc -seqfile " +
           "#{$testdata}gt_cds_descrange_multi_fail_1.fas " +
           "#{$testdata}gt_cds_test_descrange.in", :retval => 1)
  grep last_stderr, "does contain multiple sequences with ID"
end

Name "gt cds test (multi description fail 2)"
Keywords "gt_cds usedesc"
Test do
  run_test("#{$bin}gt cds -usedesc -seqfile " +
           "#{$testdata}gt_cds_descrange_multi_fail_2.fas " +
           "#{$testdata}gt_cds_test_descrange.in", :retval => 1)
  grep last_stderr, "does contain multiple sequences with ID"
end

Name "gt cds test (wrong ID)"
Keywords "gt_cds usedesc"
Test do
  run_test("#{$bin}gt cds -usedesc -seqfile " +
           "#{$testdata}gt_cds_descrange_wrong_id.fas " +
           "#{$testdata}gt_cds_test_descrange.in", :retval => 1)
  grep last_stderr, "does not contain a sequence with ID"
end

Name "gt cds test (wrong range)"
Keywords "gt_cds usedesc"
Test do
  run_test("#{$bin}gt cds -usedesc -seqfile " +
           "#{$testdata}gt_cds_descrange_wrong_range.fas " +
           "#{$testdata}gt_cds_test_descrange.in", :retval => 1)
  grep last_stderr, "cannot find sequence ID"
end

Name "gt cds test (-startcodon no -finalstopcodon no)"
Keywords "gt_cds"
Test do
  run_test "#{$bin}gt cds -startcodon no -finalstopcodon no -seqfile " +
           "#{$testdata}U89959_genomic.fas " +
           "#{$testdata}gt_cds_nostartcodon_nofinalstopcodon.in"
  run "diff #{last_stdout} " +
      "#{$testdata}gt_cds_nostartcodon_nofinalstopcodon.out"
end

Name "gt cds test (nGASP)"
Keywords "gt_cds nGASP"
Test do
  run_test "#{$bin}gt cds -startcodon yes -finalstopcodon no -minorflen 64 " +
           "-seqfile #{$testdata}nGASP/III.fas -usedesc " +
           "#{$testdata}nGASP/resIII.gff3"
  run "diff #{last_stdout} #{$testdata}nGASP/resIIIcds.gff3"
end

Name "gt cds test (U89959)"
Keywords "gt_cds"
Test do
  run_test "#{$bin}gt cds -seqfile #{$testdata}U89959_genomic.fas " +
           "#{$testdata}U89959_csas.gff3"
  run      "diff #{last_stdout} #{$testdata}U89959_cds.gff3"
end

Name "gt cds test (not sorted)"
Keywords "gt_cds"
Test do
  run_test "#{$bin}gt cds -seqfile #{$testdata}U89959_genomic.fas " +
           "#{$testdata}not_sorted.gff3", :retval => 1
  grep last_stderr, "is not sorted"
end

if $gttestdata then
  Name "gt cds bug"
  Keywords "gt_cds"
  Test do
    run_test "#{$bin}gt cds -startcodon yes -minorflen 1 " +
             "-seqfile #{$gttestdata}cds/marker_region.fas " +
             "#{$gttestdata}cds/marker_bug.gff3"
    run "diff #{last_stdout} #{$gttestdata}cds/marker_bug.out"
  end
end
