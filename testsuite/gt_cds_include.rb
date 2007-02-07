1.upto(14) do |i|
  Name "gt cds test #{i}"
  Keywords "gt_cds"
  Test do
    run_test "#{$bin}gt cds #{$testdata}gt_cds_test_#{i}.in #{$testdata}gt_cds_test_#{i}.fas"
    run "diff #{$last_stdout} #{$testdata}/gt_cds_test_#{i}.out"
  end
end
