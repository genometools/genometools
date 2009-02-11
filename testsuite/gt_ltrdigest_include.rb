if $gttestdata then
  Name "gt ltrdigest missing input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest unsorted input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3 #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest missing input FASTA"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest corrupt input FASTA"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/corrupt_input.fas", :retval => 1
  end

  Name "gt ltrdigest missing tRNA library but -trna given"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -trna #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest corrupt tRNA library"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -trna #{$gttestdata}/ltrdigest/corrupt_trna.fas #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest corrupt pHMM"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/corrupt.hmm -- #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest HMM list not properly closed (--)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/corrupt.hmm #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest tRNA implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pbsalilen 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pbsoffset 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pbstrnaoffset 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest pHMM implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pdomevalcutoff 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -pdomevalcutoff 2.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
  end

  Name "gt ltrdigest PPT HMM parameters (background distribution)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -ppttprob 0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -ppttprob 1.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -ppttprob -0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -ppttprob 0.1 -pptaprob 0.1 -pptgprob 0.2 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -ppttprob 0.25 -pptaprob 0.25 -pptgprob 0.25 -pptcprob 0.25 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 0
    run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (U-box U frequency)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pptuprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptuprob 0.0 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 0
    run_test "#{$bin}gt ltrdigest -pptuprob 0.91 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 0
    run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (PPT R/Y frequencies)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pptrprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptyprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptrprob 0.97 -pptyprob 0.03 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 0
    run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
    run_test "#{$bin}gt ltrdigest -pptrprob 0.6 -pptyprob 0.4 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz" +
             " ", :retval => 0
    run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3", :retval => 1
  end

  # positive test for all D.mel chromosomes -> must match reference
  chromosomes_dmel = ["2L","2R","3L","3R","4","X"]
  chromosomes_dmel.each do |chr|
  if $arguments["hmmer"] then
    Name "gt ltrdigest D. melanogaster chromosome #{chr} basic test w/ RT"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -threads 2 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}/ltrdigest/hmms/RVT_1_fs.hmm --  #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0, :maxtime => 12000
      run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref.gff3"
    end
  else
    Name "gt ltrdigest D. mel. chromosome #{chr} basic test, no HMM"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz",\
     :retval => 0, :maxtime => 500
      run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref_noHMM.gff3"
    end
  end
end

# XXX:disabled for now due to unexplained memleak
#  Name "gt ltrdigest GFF and FASTA do not match"
#  Keywords "gt_ltrdigest"
#  Test do
#    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
#    grep($last_stderr, "exceeds sequence boundaries!")
#  end
end
