Name "sequence buffer: EMBL 2-character code missing"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test1.embl", :retval => 1
  grep(last_stderr, "2-character line code not found in line 2")
end

Name "sequence buffer: EMBL file type unguessable"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test2.embl", :retval => 1
  grep(last_stderr, "cannot guess file type of file")
end

Name "sequence buffer: EMBL missing blanks"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test3.embl", :retval => 1
  grep(last_stderr, "3 blanks expected between")
end

Name "sequence buffer: EMBL unterminated sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test4.embl", :retval => 1
  grep(last_stderr, "unterminated sequence")
end

Name "sequence buffer: EMBL multi-line description"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test5.embl"
  run "diff #{last_stdout} #{$testdata}embl_test5.fas"
end

Name "sequence buffer: EMBL empty sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}embl_test6.embl", :retval => 1
  grep(last_stderr, "sequence 0 is empty")
end

Name "sequence buffer: EMBL files mixed with GenBank"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}Random.embl " + \
           "#{$testdata}Random.gbk", :retval => 1
end

Name "sequence buffer: EMBL files mixed with FASTA"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}Random.embl " + \
           "#{$testdata}Random.fna", :retval => 1
end

allfiles = ["Atinsert",
            "Duplicate",
            "Random-Small",
            "Random",
            "Random159",
            "Random160",
            "RandomN",
            "TTT-small",
            "trna_glutamine"]

allfiles.each do |file|
  Name "sequence buffer: check EMBL <-> FASTA #{file}"
  Keywords "gt_convertseq sequencebuffer"
  Test do
    run_test "#{$bin}gt convertseq #{$testdata}#{file}.fna " + \
             "| grep -v '>'  > #{file}_out1"
    run_test "#{$bin}gt convertseq #{$testdata}#{file}.embl " + \
             "| grep -v '>' > #{file}_out2"
    run "diff -i #{file}_out1 #{file}_out2"
  end

  Name "sequence buffer: check EMBL effectivelength #{file}"
  Keywords "gt_convertseq sequencebuffer filelengthvalues"
  Test do
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.embl"
    efflength1 = File.open(last_stderr).read.match(/ \d+\/(\d+)$/)[1].to_i
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.fna"
    efflength2 = File.open(last_stderr).read.match(/ \d+\/(\d+)$/)[1].to_i
    raise if efflength1 != efflength2
  end

  Name "sequence buffer: check EMBL rawfilelength #{file}"
  Keywords "gt_convertseq sequencebuffer filelengthvalues"
  Test do
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.embl"
    rawlength = File.open(last_stderr).read.match(/ (\d+)\/\d+$/)[1].to_i
    realsize = File.size("#{$testdata}#{file}.embl")
    raise if rawlength != realsize
  end
end

Name "sequence buffer: GenBank empty sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test1.gbk", \
           :retval => 1
  grep(last_stderr, "sequence 0 is empty")
end

Name "sequence buffer: GenBank file type unguessable"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test3.gbk", \
           :retval => 1
  grep(last_stderr, "cannot guess file type of file")
end

Name "sequence buffer: GenBank unterminated sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test4.gbk", \
           :retval => 1
  grep(last_stderr, "only terminators")
end

Name "sequence buffer: GenBank missing sequence offset numbers"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test6.gbk", \
           :retval => 1
  grep(last_stderr, "sequence offset")
end

Name "sequence buffer: GenBank blank missing"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test5.gbk", \
           :retval => 1
  grep(last_stderr, "blank expected")
end

Name "sequence buffer: GenBank text before LOCUS"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test2.gbk"
  run "diff #{last_stdout} #{$testdata}genbank_test2.fas"
end

Name "sequence buffer: GenBank missing DEFINITION line"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test2.gbk"
end

Name "sequence buffer: GenBank multiple DEFINITION lines"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}genbank_test8.gbk", \
           :retval => 1
  grep(last_stderr, "encountered another DEFINITION line")
end

Name "sequence buffer: GenBank files mixed with EMBL"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}Random.gbk " + \
           "#{$testdata}Random.embl"
end

Name "sequence buffer: GenBank files mixed with FASTA"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}Random.gbk " + \
           "#{$testdata}Random.fna"
end

gbfiles = ["Atinsert",
           "Duplicate",
           "Random-Small",
           "Random",
           "Random159",
           "Random160",
           "RandomN",
           "TTT-small",
           "trna_glutamine"]

gbfiles.each do |file|
  Name "sequence buffer: check GenBank <-> FASTA #{file}"
  Keywords "gt_convertseq sequencebuffer"
  Test do
    run_test "#{$bin}gt convertseq #{$testdata}#{file}.fna " + \
             "| grep -v '>' > #{file}_out1"
    run_test "#{$bin}gt convertseq #{$testdata}#{file}.gbk " + \
             "| grep -v '>' > #{file}_out2"
    run "diff -i #{file}_out1 #{file}_out2"
  end

  Name "sequence buffer: check GenBank effectivelength #{file}"
  Keywords "gt_convertseq sequencebuffer filelengthvalues"
  Test do
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.gbk"
    efflength1 = File.open(last_stderr).read.match(/ \d+\/(\d+)$/)[1].to_i
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.fna"
    efflength2 = File.open(last_stderr).read.match(/ \d+\/(\d+)$/)[1].to_i
    raise if efflength1 != efflength2
  end

  Name "sequence buffer: check GenBank rawfilelength #{file}"
  Keywords "gt_convertseq sequencebuffer filelengthvalues"
  Test do
    run_test "#{$bin}gt convertseq -noseq -showfilelengthvalues " + \
             "#{$testdata}#{file}.gbk"
    rawlength = File.open(last_stderr).read.match(/ (\d+)\/\d+$/)[1].to_i
    realsize = File.size("#{$testdata}#{file}.gbk")
    raise if rawlength != realsize
  end
end

Name "sequence buffer: FastQ success"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}test1.fastq"
  run "diff -i #{last_stdout} #{$testdata}test1.fasta"
end

Name "sequence buffer: FastQ success, seq > buffersize"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt dev readreads -fasta #{$testdata}fastq_long.fastq > ref.fasta"
  run_test "#{$bin}gt convertseq #{$testdata}fastq_long.fastq"
  run "diff -i #{last_stdout} ref.fasta"
end

Name "sequence buffer: FastQ non-FASTQ file"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}eden.gff3", \
           :retval => 1
  grep(last_stderr, /unknown file contents/)
end

Name "sequence buffer: FastQ invalid block start"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}test2_wrong_begin.fastq", \
           :retval => 1
  grep(last_stderr, /unknown file contents/)
end

Name "sequence buffer: FastQ different seqnames"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq " + \
           "#{$testdata}test3_different_seqnames.fastq", \
           :retval => 1
  grep(last_stderr, "sequence description 'HWI-EAS306_9_FC305MP_6_1_1331" + \
                     "_1843' is not equal to qualities description 'HWI-EAS3" +\
                     "06_9_FC305MP_6_1_1331' in line")
end

Name "sequence buffer: FastQ different seqlengths 1"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq " + \
           "#{$testdata}test4_different_seqlengths.fastq", \
           :retval => 1
  grep(last_stderr, "lengths of character sequence and qualities sequence " + \
                     "differ")
end

Name "sequence buffer: FastQ different seqlengths 2"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq " + \
           "#{$testdata}test9_uneven_length.fastq", \
           :retval => 1
  grep(last_stderr, "qualities string of sequence length 33 is not ended " + \
                     "by newline")
end

Name "sequence buffer: FastQ tricky"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt dev readreads -fasta " + \
           "#{$testdata}test5_tricky.fastq > ref.fas"
  run_test "#{$bin}gt convertseq " + \
           "#{$testdata}test5_tricky.fastq"
  run "diff -i #{last_stdout} ref.fas"
end

Name "sequence buffer: FastQ empty sequence"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}test7_empty_seq.fastq", \
           :retval => 1
  grep(last_stderr, /empty sequence/)
end

Name "sequence buffer: FastQ premature end"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}test6_premature_end.fastq", \
           :retval => 1
  grep(last_stderr, /premature end/)
end

Name "sequence buffer: FastQ multiline"
Keywords "gt_convertseq sequencebuffer fastq"
Test do
  run_test "#{$bin}gt convertseq #{$testdata}test10_multiline.fastq"
end
