Name "gt gtf_to_gff3 test"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test.gtf"
  run "diff #{last_stdout} #{$testdata}gt_gtf_to_gff3_test.gff3"
end

Name "gt gtf_to_gff3 test (lone stop codon)"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test_fail_no_flank.gtf", :retval => 1
  grep(last_stderr, /with no flanking CDS/)
end

Name "gt gtf_to_gff3 test (inconsistent strand, gene level)"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test_inconsistent_strand2.gtf", :retval => 1
  grep(last_stderr, /transcript on strand \+ encountered, but the parent gene .OR4F16. has strand/)
end

Name "gt gtf_to_gff3 test (inconsistent strand, transcript level)"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test_inconsistent_strand1.gtf", :retval => 1
  grep(last_stderr, /feature "NR_024540" on line 5 has strand \+, but the parent transcript has strand \-/)
end

Name "gt gtf_to_gff3 test (stop codon included in CDS)"
Keywords "gt_gtf_to_gff3"
Test do
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds.gtf", :retval => 1
  grep(last_stderr, /is contained in CDS in line 7/)
  run_test "#{$bin}gt gtf_to_gff3 -tidy #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds.gtf"
  grep(last_stderr, /is contained in CDS in line 7/)
  run("diff #{last_stdout} #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds.gff3")
  run_test "#{$bin}gt gtf_to_gff3 #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds2.gtf", :retval => 1
  grep(last_stderr, /is contained in CDS in line 21/)
  run_test "#{$bin}gt gtf_to_gff3 -tidy #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds2.gtf"
  grep(last_stderr, /is contained in CDS in line 21/)
  run("diff #{last_stdout} #{$testdata}gt_gtf_to_gff3_test_stop_codon_in_cds2.gff3")
end

if $gttestdata then
  Name "gt gtf_to_gff3 test D. melanogaster"
  Keywords "gt_gtf_to_gff3 large_gtf"
  Test do
    run "#{$bin}gt gff3 -sort #{$gttestdata}gff3/Drosophila_melanogaster.BDGP5.4.50.gff3 " +
        "> ref_sorted.gff3"
    run_test("#{$bin}gt gtf_to_gff3 " +
             "#{$gttestdata}gtf/Drosophila_melanogaster.BDGP5.4.50.gtf " +
             "| #{$bin}gt gff3 -sort ", :maxtime => 3600)
    run "diff #{last_stdout} " +
        "ref_sorted.gff3"
  end

  Name "gt gtf_to_gff3 test (GENCODE)"
  Keywords "gt_gtf_to_gff3 large_gtf"
  Test do
    run_test "#{$bin}gt gtf_to_gff3 -tidy " +
             "#{$gttestdata}gtf/gencode.v11.annotation.gtf.gz " +
             "| #{$bin}gt gff3 -sort -tidy", :maxtime => 3600
    run "diff #{last_stdout} #{$gttestdata}gff3/gencode.v11.annotation.gff3"
  end
end
