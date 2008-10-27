if $gttestdata then
  Name "gt ltrdigest corrupt input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz < thisisnotgff", :retval => 1
  end

  Name "gt ltrdigest unsorted input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest missing input FASTA"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest corrupt input FASTA"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -seqfile #{$gttestdata}/ltrdigest/corrupt_input.fas #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest missing tRNA library but -trna given"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -trna -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted ", :retval => 1
  end

  Name "gt ltrdigest corrupt tRNA library"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -trna #{$gttestdata}/ltrdigest/corrupt_trna.fas -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest corrupt pHMM"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/corrupt.hmm -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest HMM list not properly closed (--)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/corrupt.hmm -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted ", :retval => 1
  end

  Name "gt ltrdigest tRNA implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pbsalilen 10 20 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
    run_test "#{$bin}gt ltrdigest -pbsoffset 10 20 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
    run_test "#{$bin}gt ltrdigest -pbstrnaoffset 10 20 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted", :retval => 1
  end

  Name "gt ltrdigest pHMM implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -pdomevalcutoff 0.2 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted " +
             " ", :retval => 1
    run_test "#{$bin}gt ltrdigest -pdomevalcutoff 2.2 -seqfile #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted" +
             " ", :retval => 1
  end

  Name "gt ltrdigest GFF and FASTA do not match"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest -seqfile #{$gttestdata}/ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted ", :retval => 1
    grep($last_stderr, "exceeds sequence boundaries!")
  end

  Name "gt ltrdigest D. mel. all chromosomes, mapping failure"
  Keywords "gt_ltrdigest"
  Test do
      txt = "function mapping(sequence_region)
               if sequence_region == '2L' then
                 return 'foo'
               end
               return \"#{$gttestdata}/ltrdigest/\"..sequence_region..\"_genomic_dmel_RELEASE3-1.FASTA.gz\"
             end"
      File.open("mapping.lua","w+") do |file|
        file.write(txt)
      end
      run_test "#{$bin}gt ltrdigest -regionmapping mapping.lua -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_all_chromosomes.sorted.gff3 ",\
               :retval => 1
      grep($last_stderr, " does not exist or is not readable")
  end

  Name "gt ltrdigest D. mel. all chromosomes, with mapping, no HMM"
  Keywords "gt_ltrdigest"
  Test do
    txt = "function mapping(sequence_region)
             return \"#{$gttestdata}/ltrdigest/\"..sequence_region..\"_genomic_dmel_RELEASE3-1.FASTA.gz\"
           end"
    File.open("mapping.lua","w+") do |file|
      file.write(txt)
    end
    run_test "#{$bin}gt ltrdigest -regionmapping mapping.lua -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_all_chromosomes.sorted.gff3 ", :retval => 0, :maxtime => 4500
    run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/dmel_all_chromosomes.ref.gff3"
  end

  # positive test for all D.mel chromosomes -> must match reference
  chromosomes_dmel = ["2L","2R","3L","3R","4","X"]
  chromosomes_dmel.each do |chr|
    Name "gt ltrdigest D. mel. chromosome single region #{chr} w/ mapping"
    Keywords "gt_ltrdigest"
    Test do
      txt = "function mapping(sequence_region)
               if sequence_region == 'seq0' then
                 return \"#{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz\"
               else
                 return ''
               end
             end"
      File.open("mapping.lua","w+") do |file|
        file.write(txt)
      end
      run_test "#{$bin}gt ltrdigest -outfileprefix out#{chr} -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa -regionmapping mapping.lua #{$gttestdata}/ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted", :retval => 0, :maxtime => 1000
      run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref_noHMM.gff3"
      run "test -e out#{chr}_conditions.csv"
    end
  end

  chromosomes_dmel.each do |chr|
    if $arguments["hmmer"] then
      Name "gt ltrdigest D. melanogaster chromosome #{chr} basic test w/ RT"
      Keywords "gt_ltrdigest"
      Test do
        run_test "#{$bin}gt ltrdigest -threads 2 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}/ltrdigest/hmms/RVT_1_fs.hmm -seqfile #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz  #{$gttestdata}/ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted ",\
         :retval => 0, :maxtime => 2000
        run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref.gff3"
      end
    else
      Name "gt ltrdigest D. mel. chromosome #{chr} basic test, no HMM"
      Keywords "gt_ltrdigest"
      Test do
          run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa -seqfile #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}/ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted ",\
         :retval => 0, :maxtime => 500
          run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref_noHMM.gff3"
      end
    end
  end
end
