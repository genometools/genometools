1.upto(14) do |i|
  Name "gt cds test #{i}"
  Keywords "gt_cds"
  Test do
    run_test "#{$bin}gt cds -seqfile #{$testdata}gt_cds_test_#{i}.fas #{$testdata}gt_cds_test_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}/gt_cds_test_#{i}.out"
  end
end

1.upto(14) do |i|
  Name "gt cds test #{i} (-usedesc)"
  Keywords "gt_cds usedesc"
  Test do
    run_test "#{$bin}gt cds -usedesc -seqfile #{$testdata}gt_cds_test_#{i}.fas #{$testdata}gt_cds_test_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}/gt_cds_test_#{i}.out"
  end
end

Name "gt cds test (description range)"
Keywords "gt_cds usedesc"
Test do
  run_test "#{$bin}gt cds -usedesc -seqfile " +
           "#{$testdata}gt_cds_test_descrange.fas " +
           "#{$testdata}gt_cds_test_descrange.in"
  run "diff #{$last_stdout} #{$testdata}/gt_cds_test_descrange.out"
end

if $gttestdata then
  Name "gt cds bug"
  Keywords "gt_cds"
  Test do
    run_test "#{$bin}gt cds -seqfile #{$gttestdata}cds/marker_region.fas " +
             "#{$gttestdata}/cds/marker_bug.gff3"
    run "diff #{$last_stdout} #{$gttestdata}/cds/marker_bug.out"
  end
end
