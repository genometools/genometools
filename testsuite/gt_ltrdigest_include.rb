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

  # positive test for all D.mel chromosomes -> must match reference
  chromosomes_dmel = ["2L","2R","3L","3R","4","X"]
  chromosomes_dmel.each do |chr|
  if $arguments["hmmer"] then
    Name "gt ltrdigest D. melanogaster chromosome #{chr} basic test w/ RT"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}/ltrdigest/hmms/RVT_1_fs.hmm --  #{$gttestdata}/ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz",\
       :retval => 0, :maxtime => 500
      run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref.gff3"
    end
  else
    Name "gt ltrdigest D. mel. chromosome #{chr} basic test, no HMM"
    Keywords "gt_ltrdigest"
    Test do
        run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz",\
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
