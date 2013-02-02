# simple reverse complement
def revcomp(seq)
  seq.tr("aAcCgGtTnN","tTgGcCaSAnN").reverse
end

# simple DNA->AA translation
def translate_dna(sequence)
  genetic_code = {
  'ATA'=>'I', 'ATC'=>'I', 'ATT'=>'I', 'ATG'=>'M',
  'ACA'=>'T', 'ACC'=>'T', 'ACG'=>'T', 'ACT'=>'T',
  'AAC'=>'N', 'AAT'=>'N', 'AAA'=>'K', 'AAG'=>'K',
  'AGC'=>'S', 'AGT'=>'S', 'AGA'=>'R', 'AGG'=>'R',
  'CTA'=>'L', 'CTC'=>'L', 'CTG'=>'L', 'CTT'=>'L',
  'CCA'=>'P', 'CCC'=>'P', 'CCG'=>'P', 'CCT'=>'P',
  'CAC'=>'H', 'CAT'=>'H', 'CAA'=>'Q', 'CAG'=>'Q',
  'CGA'=>'R', 'CGC'=>'R', 'CGG'=>'R', 'CGT'=>'R',
  'GTA'=>'V', 'GTC'=>'V', 'GTG'=>'V', 'GTT'=>'V',
  'GCA'=>'A', 'GCC'=>'A', 'GCG'=>'A', 'GCT'=>'A',
  'GAC'=>'D', 'GAT'=>'D', 'GAA'=>'E', 'GAG'=>'E',
  'GGA'=>'G', 'GGC'=>'G', 'GGG'=>'G', 'GGT'=>'G',
  'TCA'=>'S', 'TCC'=>'S', 'TCG'=>'S', 'TCT'=>'S',
  'TTC'=>'F', 'TTT'=>'F', 'TTA'=>'L', 'TTG'=>'L',
  'TAC'=>'Y', 'TAT'=>'Y', 'TAA'=>'X', 'TAG'=>'X',
  'TGC'=>'C', 'TGT'=>'C', 'TGA'=>'X', 'TGG'=>'W',
  }
  proteinseq = ""
  (0..sequence.length).step(3) do |i|
    codon = sequence[i..i+2].upcase
    next if codon.length != 3
    if genetic_code.has_key?(codon) then
      proteinseq += genetic_code[codon]
    else
      raise "unknown codon '#{codon}'"
    end
  end
  proteinseq
end

# returns a hash indexed by element coordinates
# containing an array of DNA sequences (one per hit fragment)
def get_protdom_dna_seqs(gff3filename, chr)
  seqs = Hash.new
  File.open(gff3filename) do |file|
    mainstart, mainstop = 0, 0
    mainstrand = ''
    file.readlines.each do |ln|
      _, _, type, start, stop, _, strand = ln.split("\t")
      if type == "LTR_retrotransposon" then
        mainstart = start.to_i
        mainstop = stop.to_i
        mainstrand = strand
      elsif type == "protein_match" then
        run_test "#{$bin}gt extractseq -frompos " + \
                 "#{start} -topos " + \
                 "#{stop} " + \
                 "#{$gttestdata}ltrdigest/#{chr}_genomic_dmel_" + \
                 "RELEASE3-1.FASTA.gz > tmp.fas"
        outa = File.open("tmp.fas").read.split("\n")
        outa.shift
        if seqs[[mainstart, mainstop]].nil? then
          seqs[[mainstart, mainstop]] = []
        end
        out = outa.collect{|l| l.chomp.downcase}.join
        # we want to provide DNA sequences in reading direction
        if mainstrand == '-' then
          seqs[[mainstart, mainstop]].insert(0,(revcomp(out)))
        else
          seqs[[mainstart, mainstop]].push(out)
        end
      end
    end
  end
  seqs
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
          run_test "#{$bin}gt extractseq -frompos " + \
                   "#{coords[feat][[start,stop]][0]} -topos " + \
                   "#{coords[feat][[start,stop]][1]} " + \
                   "#{$gttestdata}ltrdigest/#{chr}_genomic_dmel_" + \
                   "RELEASE3-1.FASTA.gz > tmp.fas"
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

def load_seqfile(filename)
  contents = {}
  File.open(filename) do |file|
    (farr = file.read.split(">")).each do |fitem|
      seqa = fitem.split("\n")
      fline = seqa.shift
      next if fline.nil?
      if m = /_(\d+)_(\d+)$/.match(fline) then
        start, stop = m[1].to_i, m[2].to_i
        seq = seqa.collect{|l| l.chomp.downcase}.join
        contents[[start,stop]] = seq
      end
    end
  end
  contents
end

def check_amino_acid_output(gff3file, chr, pdoms)
  pdoms.each do |dom|
    dnaseqfile = load_seqfile("result#{chr}_pdom_#{dom}.fas")
    aaseqfile  = load_seqfile("result#{chr}_pdom_#{dom}_aa.fas")
    seqs = get_protdom_dna_seqs(gff3file, chr)
    if dnaseqfile.keys.length != aaseqfile.keys.length then
      raise "number of DNA and AA sequences is not equal"
    end
    if (dnaseqfile.keys - aaseqfile.keys).length != 0 then
      raise "DNA and AA sequence keys are not identical"
    end
    dnaseqfile.keys.each do |key|
      if seqs[key].nil? or seqs[key].length == 0 then
        raise "no protein domain sequences extracted from GFF3+sequence " + \
              "for element #{key.inspect} of chromosome #{chr}"
      elsif seqs[key].length == 1 then
        # hit was not fragmented, we can check directly
        # first, DNA sequence in the output file must match the GFF coordinates
        unless seqs[key][0] == dnaseqfile[key] then
          raise "DNA sequence from FASTA file \n" + \
                "'#{seqs[key][0]}' does " + \
                "not match the sequence referred to in the GFF3 file \n" + \
                "'#{dnaseqfile[key]}' in element " + \
                "#{key.inspect}"
        end
        # also, translated DNA seq must match the AA sequence from FASTA out
        if translate_dna(seqs[key][0]) != aaseqfile[key].upcase then
          raise "translated sequence \n'#{translate_dna(seqs[key][0])}' " + \
                "is not equal to sequence \n'#{aaseqfile[key].upcase}' in " + \
                "element #{key.inspect}"
        end
      else
        # number of fragments > 1
        # translate each fragment's DNA sequence separately and concatenate
        # into one => this must match the AA sequence output by LTRdigest
        translated_seqs = seqs[key].collect{|s| translate_dna(s)}.join
        if translated_seqs != aaseqfile[key].upcase then
          raise "sequence translated from #{seqs[key].length} frags\n" + \
                "'#{translated_seqs}' " + \
                "is not equal to sequence \n'#{aaseqfile[key].upcase}' in " + \
                "element #{key.inspect}"
        end
      end
    end
  end
end

if $gttestdata then
  Name "gt ltrdigest missing input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt ltrdigest #{$gttestdata}ltrdigest/2L_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
  end

  Name "gt ltrdigest unsorted input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /is not sorted/)
  end

  Name "gt ltrdigest wrong sequence regions in input GFF"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}ltrdigest/dmel_test_Run9_4_wrong_seqid.gff3 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, "Feature 'LTR_retrotransposon1' on line 6 has invalid region identifier,must be 'seqX' with X being a sequence number, but was 'chr4'")
  end

  Name "gt ltrdigest missing input sequence"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted", :retval => 1
    grep(last_stderr, /missing argument/)
  end

  Name "gt ltrdigest missing tRNA library but -trna given"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -trnas -outfileprefix foo #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /missing argument to option/)
  end

  Name "gt ltrdigest corrupt tRNA library"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}ltrdigest/corrupt_trna.fas #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /cannot guess file type of file/)
  end

  Name "gt ltrdigest tRNA implied options"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pbsalilen 10 20 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /option "-pbsalilen" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -pbsoffset 10 20 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /option "-pbsoffset" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -pbstrnaoffset 10 20 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /option "-pbstrnaoffset" requires option "-trnas"/)
    run_test "#{$bin}gt ltrdigest -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{last_stdout} #{$gttestdata}ltrdigest/4_ref_noHMM.gff3"
  end

  if $arguments["hmmer"] then
    Name "gt ltrdigest corrupt pHMM"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}ltrdigest/corrupt.hmm -- #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep(last_stderr, /Failed to open HMM file/)
    end

    Name "gt ltrdigest HMM list not properly closed (--)"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep(last_stderr, /missing argument/)
    end

    Name "gt ltrdigest pHMM implied options"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -pdomevalcutoff 0.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep(last_stderr, /option "-pdomevalcutoff" requires option "-hmms"/)
      run_test "#{$bin}gt ltrdigest -pdomevalcutoff 2.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
      grep(last_stderr, /argument to option "-pdomevalcutoff" must be a floating point value <= 1.000000/)
    end

    Name "gt ltrdigest use of deprecated '-threads' switch"
    Keywords "gt_ltrdigest"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
      run_test "#{$bin}gt ltrdigest -threads 2 -outfileprefix result4 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm --  #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0, :maxtime => 12000
      grep(last_stderr, /option is deprecated. Please use/)
    end
  end

  Name "gt ltrdigest PPT HMM parameters (background distribution)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    # probabilities do not add up to 1, sum > 0
    run_test "#{$bin}gt ltrdigest -ppttprob 0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # probability > 1
    run_test "#{$bin}gt ltrdigest -ppttprob 1.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # negative probability
    run_test "#{$bin}gt ltrdigest -ppttprob -0.3 -pptaprob 0.3 -pptgprob 0.9 -pptcprob 0.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # probabilities do not add up to 1, sum < 0
    run_test "#{$bin}gt ltrdigest -ppttprob 0.1 -pptaprob 0.1 -pptgprob 0.2 -pptcprob 0.2 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    # positive test
    run_test "#{$bin}gt ltrdigest -ppttprob 0.25 -pptaprob 0.25 -pptgprob 0.25 -pptcprob 0.25 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{last_stdout} #{$gttestdata}ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (U-box U frequency)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pptuprob 1.3 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptuprob 0.0 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    run_test "#{$bin}gt ltrdigest -pptuprob 0.91 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{last_stdout} #{$gttestdata}ltrdigest/4_ref_noHMM.gff3"
  end

  Name "gt ltrdigest PPT HMM parameters (PPT R/Y frequencies)"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest -pptrprob 1.3 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptyprob 1.3 #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    run_test "#{$bin}gt ltrdigest -pptrprob 0.97 -pptyprob 0.03 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{last_stdout} #{$gttestdata}ltrdigest/4_ref_noHMM.gff3"
    run_test "#{$bin}gt ltrdigest -pptrprob 0.6 -pptyprob 0.4 -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0
    #run "diff #{last_stdout} #{$gttestdata}ltrdigest/4_ref_noHMM.gff3", :retval => 1
  end

  Name "gt ltrdigest GFF and sequence do not match"
  Keywords "gt_ltrdigest"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt ltrdigest #{$gttestdata}ltrdigest/dmel_test_Run9_2L.gff3.sorted 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, "exceeds sequence boundaries!")
  end

  # positive test for all D.mel chromosomes -> must match reference
  chromosomes_dmel = ["2L","2R","3L","3R","4","X"]
  chromosomes_dmel.each do |chr|
    if $arguments["hmmer"] then
      Name "gt ltrdigest D. melanogaster chromosome #{chr} basic test w/ RT"
      Keywords "gt_ltrdigest"
      Test do
        run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", :maxtime => 600
        run_test "#{$bin}gt ltrdigest -outfileprefix result#{chr} -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa -hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm --  #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 0, :maxtime => 12000
        check_ppt_pbs(last_stdout, chr)
        #run "diff #{last_stdout} #{$gttestdata}ltrdigest/#{chr}_ref.gff3"
      end
    else
      Name "gt ltrdigest D. mel. chromosome #{chr} basic test, no HMM"
      Keywords "gt_ltrdigest"
      Test do
        run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", :maxtime => 600
        run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result#{chr} -trnas #{$gttestdata}ltrdigest/Dm-tRNAs-uniq.fa #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz",\
       :retval => 0, :maxtime => 700
        check_ppt_pbs(last_stdout, chr)
        #run "diff #{last_stdout} #{$gttestdata}ltrdigest/#{chr}_ref_noHMM.gff3"
      end
    end
  end

  chromosomes_dmel.each do |chr|
    if $arguments["hmmer"] then
      Name "gt ltrdigest D. melanogaster chromosome #{chr} AAseq out"
      Keywords "gt_ltrdigest aminoacidout"
      Test do
        run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v " + \
                 "-db #{$gttestdata}ltrdigest/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", \
                 :maxtime => 600
        run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result#{chr} " + \
                 "-hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm " + \
                 "-aaout yes " + \
                 " #{$gttestdata}ltrdigest/dmel_test_Run9_#{chr}.gff3.sorted " + \
                 " #{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", \
                 :retval => 0, :maxtime => 12000
        check_amino_acid_output(last_stdout, chr, ["RVT_1"])
      end
    end
  end

  if $arguments["hmmer"] then
    Name "gt ltrdigest -aliout"
    Keywords "gt_ltrdigest aminoacidout aliout"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v " + \
               "-db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :maxtime => 600
      run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result4 " + \
               "-hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm -- " + \
               "#{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted " + \
               "4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :retval => 0, :maxtime => 12000
      if File.exists?("result4_pdom_RVT_1.ali") then
        raise TestFailed, "file \"result4_pdom_RVT_1.ali\" should not exist"
      end
      run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result4 " + \
               "-hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm " + \
               "-aliout yes " + \
               "#{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted " + \
               "4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :retval => 0, :maxtime => 12000
      if !File.exists?("result4_pdom_RVT_1.ali") then
        raise TestFailed, "file \"result4_pdom_RVT_1.ali\" does not exist"
      end
    end

    Name "gt ltrdigest -aaout"
    Keywords "gt_ltrdigest aminoacidout aaout"
    Test do
      run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v " + \
               "-db #{$gttestdata}ltrdigest/4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :maxtime => 600
      run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result4 " + \
               "-hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm -- " + \
               "#{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted " + \
               "4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :retval => 0, :maxtime => 12000
      if File.exists?("result4_pdom_RVT_1_aa.fas") then
        raise TestFailed, "file \"result4_pdom_RVT_1_aa.fas\" should not exist"
      end
      run_test "#{$bin}gt -j 2 ltrdigest -outfileprefix result4 " + \
               "-hmms #{$gttestdata}ltrdigest/hmms/RVT_1.hmm " + \
               "-aaout yes " + \
               "#{$gttestdata}ltrdigest/dmel_test_Run9_4.gff3.sorted " + \
               "4_genomic_dmel_RELEASE3-1.FASTA.gz", \
               :retval => 0, :maxtime => 12000
      if !File.exists?("result4_pdom_RVT_1_aa.fas") then
        raise TestFailed, "file \"result4_pdom_RVT_1_aa.fas\" does not exist"
      end
    end
  end
end
