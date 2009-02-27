# simple reverse complement
def revcomp(seq)
  seq.tr("aAcCgGtTnN","tTgGcCaSAnN").reverse
end

# checks whether the sequences output by the LTRFileOutStream match those
# referenced in the GFF3 file
def check_ppt_pbs(gff3filename, chr)
  coords = {:ppt => {}, :pbs => {}}
  File.open(gff3filename) do |file|
    mainstart, mainstop = 0, 0
    file.readlines.each do |ln|
      _, _, type, start, stop = ln.split("\t")
      if type == "LTR_retrotransposon" then
        mainstart = start.to_i
        mainstop = stop.to_i
      elsif type == "RR_tract" then
        coords[:ppt][[mainstart,mainstop]] = [start.to_i, stop.to_i]
      elsif type == "primer_binding_site" then
        coords[:pbs][[mainstart,mainstop]] = [start.to_i, stop.to_i]
      end
    end
  end
  [:ppt,:pbs].each do |feat|
    File.open("result#{chr}_#{feat}.fas") do |file|
      (farr = file.read.split(">")).each do |fitem|
        seqa = fitem.split("\n")
        fline = seqa.shift
        next if fline.nil?
        if m = /_(\d+)_(\d+)$/.match(fline) then
          start, stop = m[1].to_i, m[2].to_i
          seq = seqa.collect{|l| l.chomp.downcase}.join
          run_test "#{$bin}gt extractseq -frompos #{coords[feat][[start,stop]][0]} -topos #{coords[feat][[start,stop]][1]} #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz > tmp.fas"
          outa = File.open("tmp.fas").read.split("\n")
          outa.shift
          out = outa.collect{|l| l.chomp.downcase}.join
          unless out == seq or revcomp(out) == seq
            raise "sequence output and GFF3 coordinates do not match " + \
                  "(#{seq} != #{out}) for the #{feat} with coords " + \
                  "#{coords[feat][[start,stop]][0]}-" + \
                  "#{coords[feat][[start,stop]][1]}, neither does reverse " + \
                  "complement!"
          end
        end
      end
    end
  end
end

if $gttestdata then
  Name "gt ltrdigest missing input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest unsorted input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /is not sorted/)
  end

  Name "gt ltrdigest wrong sequence regions in input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_4_wrong_seqid.gff3 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, "Feature 'LTR_retrotransposon1' on line 6 has invalid region identifier,must be 'seqX' with X being a sequence number, but was 'chr4'")
  end

  Name "gt ltrdigest missing input sequence"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted", :retval => 1
    grep($last_stderr, /missing argument/)
  end

  Name "gt ltrdigest missing tRNA library but -trna given"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -trnas -outfileprefix foo #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /superfluous argument/)
  end

  Name "gt ltrdigest corrupt tRNA library"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/corrupt_trna.fas #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /the first character of fasta file/)
  end

  Name "gt ltrdigest tRNA implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pbsalilen 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /option "-pbsalilen" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -pbsoffset 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /option "-pbsoffset" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -pbstrnaoffset 10 20 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, /option "-pbstrnaoffset" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
  end

  if $arguments["hmmer"] then
    Name "gt ltrdigest corrupt pHMM"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/corrupt.hmm -- #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep($last_stderr, /Failed to open HMM file/)
    end

    Name "gt ltrdigest HMM list not properly closed (--)"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}/ltrdigest/hmms/RVT_1_fs.hmm #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep($last_stderr, /missing argument/)
    end
  
    Name "gt ltrdigest pHMM implied options"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -pdomevalcutoff 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep($last_stderr, /option "-pdomevalcutoff" requires option "-hmms"/)
      run_test "#{$bin}gt ltrdigest -pdomevalcutoff 2.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep($last_stderr, /argument to option "-pdomevalcutoff" must be a floating point value <= 1.000000/)
    end
  end

  Name "gt ltrdigest PPT HMM parameters (background distribution)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    # probabilities do not add up to 1, sum > 0
    run_test "#{$bin}gt ltrdigest -ppttprob 0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # probability > 1
    run_test "#{$bin}gt ltrdigest -ppttprob 1.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # negative probability
    run_test "#{$bin}gt ltrdigest -ppttprob -0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # probabilities do not add up to 1, sum < 0
    run_test "#{$bin}gt ltrdigest -ppttprob 0.1 -pptaprob 0.1 -pptgprob 0.2 -pptcprob 0.2 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # positive test
    run_test "#{$bin}gt ltrdigest -ppttprob 0.25 -pptaprob 0.25 -pptgprob 0.25 -pptcprob 0.25 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (U-box U frequency)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pptuprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptuprob 0.0 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    run_test "#{$bin}gt ltrdigest -pptuprob 0.91 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (PPT R/Y frequencies)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pptrprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptyprob 1.3 #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptrprob 0.97 -pptyprob 0.03 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3"
    run_test "#{$bin}gt ltrdigest -pptrprob 0.6 -pptyprob 0.4 -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}/ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/4_ref_noHMM.gff3", :retval => 1
  end

  Name "gt ltrdigest GFF and sequence do not match"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}/ltrdigest/dmel_test_Run9_2L.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep($last_stderr, "exceeds sequence boundaries!")
  end
  
  # positive test for all D.mel chromosomes -> must match reference
  chromosomes_dmel = ["2L","2R","3L","3R","4","X"]
  chromosomes_dmel.each do |chr|
    if $arguments["hmmer"] then
      Name "gt ltrdigest D. melanogaster chromosome #{chr} basic test w/ RT"
      Keywords "gt_ltrdigest"
      Test do
        run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz"
        run_test "#{$bin}gt ltrdigest -outfileprefix result#{chr} -threads 2 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}/ltrdigest/hmms/RVT_1_fs.hmm --  #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0, :maxtime => 12000
        check_ppt_pbs($last_stdout, chr)       
        #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref.gff3"
      end
    else
      Name "gt ltrdigest D. mel. chromosome #{chr} basic test, no HMM"
      Keywords "gt_ltrdigest"
      Test do
        run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}/ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz"
        run_test "#{$bin}gt ltrdigest -outfileprefix result#{chr} -trnas #{$gttestdata}/ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz",\
       :retval => 0, :maxtime => 500
        check_ppt_pbs($last_stdout, chr)
        #run "diff #{$last_stdout} #{$gttestdata}/ltrdigest/#{chr}_ref_noHMM.gff3"
      end
    end
  end
end
