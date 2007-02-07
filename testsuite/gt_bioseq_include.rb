Name "gt bioseq nonexistent file"
Keywords "gt_bioseq"
Test do
  run_test("#{$bin}gt bioseq #{$testdata}nonexistent_file", :retval => 1)
end

1.upto(7) do |i|
  Name "gt bioseq fail #{i}"
  Keywords "gt_bioseq"
  Test do
    run_test("#{$bin}gt bioseq -recreate #{$testdata}gt_bioseq_fail_#{i}.fas",
             :retval => 1)
  end
end

Name "gt bioseq test 1"
Keywords "gt_bioseq"
Test do
  run_test "#{$bin}gt bioseq -recreate #{$testdata}gt_bioseq_succ_1.fas"
end

Name "gt bioseq test 2"
Keywords "gt_bioseq"
Test do
  run_test "#{$bin}gt bioseq -recreate #{$testdata}gt_bioseq_succ_2.fas"
end

Name "gt bioseq test 3"
Keywords "gt_bioseq"
Test do
  run_test "#{$bin}gt bioseq -recreate -showfasta -width 70 #{$testdata}gt_bioseq_succ_3.fas" 
  run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.fas"
end

1.upto(3) do |i|
  Name "gt bioseq test 3 out #{i}"
  Keywords "gt_bioseq"
  Test do
    run_test "#{$bin}gt bioseq -showseqnum #{i} -width 70 #{$testdata}gt_bioseq_succ_3.fas"
    run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.out#{i}"
  end
end

Name "gt bioseq test 3 out 4 fail"
Keywords "gt_bioseq"
Test do
  run_test("#{$bin}gt bioseq -showseqnum 4 #{$testdata}gt_bioseq_succ_3.fas",
           :retval => 1)
end

Name "gt bioseq test 3 stat"
Keywords "gt_bioseq"
Test do
  run_test "#{$bin}gt bioseq -stat #{$testdata}gt_bioseq_succ_3.fas"
end
