Name "gt csa test"
Keywords "gt_csa"
Test do
  run_test("#{$bin}gt csa #{$testdata}gt_csa_test_1", :retval => 1)
end

1.upto(4) do |i|
  Name "gt csa prob #{i}"
  Keywords "gt_csa"
  Test do
    run_test "#{$bin}gt csa #{$testdata}gt_csa_prob_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}gt_csa_prob_#{i}.out"
  end
end

1.upto(4) do |i|
  Name "gt -debug csa prob #{i}"
  Keywords "gt_csa"
  Test do
    run_test "#{$bin}gt -debug csa #{$testdata}gt_csa_prob_#{i}.in"
    run "tail -n +2 #{$last_stderr} | diff - #{$testdata}gt_csa_prob_#{i}.debug"
  end
end

Name "gt csa arabidopsis"
Keywords "gt_csa"
Test do
  run_test "#{$bin}gt csa #{$testdata}U89959_sas.gff3"
  run "diff #{$last_stdout} #{$testdata}U89959_csas.gff3"
end
