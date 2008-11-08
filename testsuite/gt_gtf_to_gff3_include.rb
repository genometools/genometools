Name "gt gtf_to_gff3 test"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test.gtf"
  run "diff #{$last_stdout} #{$testdata}gt_gtf_to_gff3_test.gff3"
end

if $gttestdata then
  Name "gt gtf_to_gff3 test"
  Keywords "gt_gtf_to_gff3 large_gtf"
  Test do
  run_test("#{$bin}gt gtf_to_gff3 " +
           "#{$gttestdata}gtf/Drosophila_melanogaster.BDGP5.4.50.gtf",
           :maxtime => 30)
  run("diff #{$last_stdout} " +
      "#{$gttestdata}gff3/Drosophila_melanogaster.BDGP5.4.50.gff3",
      :maxtime => 20)
  end
end
