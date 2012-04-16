rdjO = "#{$bin}gt readjoiner overlap"
rdjA = "#{$bin}gt readjoiner assembly"
spmtest = "#{$bin}gt readjoiner spmtest"
cnttest = "#{$bin}gt readjoiner cnttest"

def run_correct(input, k, c)
  run "#{$bin}gt suffixerator -mirrored -suf -lcp -ssp"+
      " -indexname reads -db #{input}"
  run "#{$bin}gt readjoiner correct -k #{k} -c #{c} -ii reads"
end

def run_prefilter(input, moreopts="")
  run "#{$bin}gt readjoiner prefilter -readset reads "+
      "-db #{input} #{moreopts}"
end

def run_overlap(minlen, moreopts="")
  run "#{$bin}gt readjoiner overlap -readset reads "+
    "-l #{minlen} #{moreopts}"
end

def run_assembly(moreopts="")
  run "#{$bin}gt readjoiner assembly -readset reads #{moreopts}"
end

def encode_reads(input, moreopts="")
  run "#{$bin}gt encseq encode -indexname reads "+
      "#{moreopts} #{input}"
end

def prepare_esa(indexname, mirrored)
  run "#{$bin}gt suffixerator -ii #{indexname} -suf -lcp -ssp "+
      "#{mirrored ? '-mirrored' : ''}"
end

Name "gt readjoiner correct"
Keywords "gt_readjoiner gt_readjoiner_correct"
Test do
  run_correct("#{$testdata}/readjoiner/errors_1.fas", 12, 2)
  run "#{$bin}gt encseq decode reads"
  run "diff #{last_stdout} #{$testdata}/readjoiner/errors_1.corrected.fas"
end

Name "gt readjoiner prefilter: correct encseq output (eqlen)"
Keywords "gt_readjoiner gt_readjoiner_prefilter"
Test do
  # first prepare contfree readset in fasta format
  run_prefilter("#{$testdata}/readjoiner/30x_800nt.fas",
                "-encseq false -fasta true -q true")
  contfree = last_stdout
  # prepare encseq using prefilter
  run_prefilter(contfree)
  run "mv reads.esq reads.prefilter.esq"
  # prepare encseq using encseq encode
  encode_reads(contfree)
  run "diff reads.esq reads.prefilter.esq"
end

Name "gt readjoiner prefilter: multiple input files (eqlen)"
Keywords "gt_readjoiner gt_readjoiner_prefilter"
Test do
  run_prefilter("#{$testdata}/readjoiner/2x3nt_1.fas "+
                "#{$testdata}/readjoiner/2x3nt_2.fas "+
                "#{$testdata}/readjoiner/2x3nt_3.fas")
end

Name "gt readjoiner prefilter: correct encseq output (not eqlen)"
Keywords "gt_readjoiner gt_readjoiner_prefilter"
Test do
  run_prefilter("#{$testdata}/readjoiner/30x_long_varlen.fas",
                "-encseq false -fasta true -q true")
  contfree = last_stdout
  run_prefilter(contfree, "")
  run "mv reads.esq reads.prefilter.esq"
  run "mv #{contfree} reads.prefiltered.fas"
  encode_reads("reads.prefiltered.fas")
  run "diff reads.esq reads.prefilter.esq"
end

Name "gt readjoiner: no crash when no spm found"
Keywords "gt_readjoiner"
Test do
   run_prefilter("#$testdata/readjoiner/tiny.fas")
   run_overlap(12)
   run "grep 'number of irreducible suffix-prefix matches = 0' #{last_stdout}"
   run_assembly
   run "grep 'no contigs' #{last_stdout}"
end

Name "gt readjoiner overlap: eqlen; minlen > readlen"
Keywords "gt_readjoiner gt_readjoiner_overlap"
Test do
  run_prefilter("#$testdata/readjoiner/tiny.fas")
  run_test "#{$bin}gt readjoiner overlap -readset reads -l 17"
end

=begin singlestrand
Name "gt readjoiner overlap: singlestrand vs mirrored mode"
Keywords "gt_readjoiner gt_readjoiner_overlap"
Test do
  run_prefilter("#$testdata/readjoiner/tiny.fas")
  run_overlap(4, "-singlestrand")
  run "#{spmtest} -test showlist -readset reads"
  spm = last_stdout
  run_test "diff #{spm} " +
    "#$testdata/readjoiner/tiny_singlestrand.spm", :retval => 0
  run_overlap(4)
  run "#{spmtest} -test showlist -readset reads"
  spm = last_stdout
  run_test "diff #{spm} " +
    "#$testdata/readjoiner/tiny_mirrored.spm", :retval => 0
end
=end

def assert_empty(fn, desc)
  fncontent = IO.read(fn)
  fncontent = fncontent.split("\n")
  # tolerate this ruby warning:
  fncontent.reject! {|line| line =~ /Insecure world writable/}
  failtest(desc + " should be empty") unless fncontent.join == ""
end

def assert_not_empty(fn, desc)
  failtest(desc + " should not be empty") unless IO.read(fn) != ""
end

Name "gt readjoiner: verbosity levels (-q/default/-v)"
Keywords "gt_readjoiner gt_readjoiner_prefilter "+
         "gt_readjoiner_overlap gt_readjoiner_assembly"
Test do
  # -- prefilter -- #

  # default mode
  run_prefilter("#$testdata/readjoiner/tiny.fas")
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # quiet mode
  run_prefilter("#$testdata/readjoiner/tiny.fas", "-q")
  assert_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # verbose mode
  run_prefilter("#$testdata/readjoiner/tiny.fas", "-v")
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")

  # -- overlap -- #

  # default mode
  run_overlap(4)
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # quiet mode
  run_overlap(4, "-q")
  assert_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # verbose mode
  run_overlap(4, "-v")
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")

  # -- assembly -- #

  # default mode
  run_assembly
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # quiet mode
  run_assembly("-q")
  assert_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
  # verbose mode
  run_assembly("-v")
  assert_not_empty(last_stdout, "stdout")
  assert_empty(last_stderr, "stderr")
end

Name "gt readjoiner: wildcards"
Keywords "gt_readjoiner gt_readjoiner_prefilter"
Test do
  run_prefilter("#$testdata/readjoiner/wildcard1.fas")
  grep last_stdout, 'reads with ambiguities = 1'
  grep last_stdout, 'number of reads in filtered readset = 1'
  run_prefilter("#$testdata/readjoiner/wildcard2.fas")
  grep last_stdout, 'reads with ambiguities = 1'
  grep last_stdout, 'number of reads in filtered readset = 1'
  run_prefilter("#$testdata/readjoiner/wildcards.fas")
  grep last_stdout, 'reads with ambiguities = 6'
  grep last_stdout, 'number of reads in filtered readset = 2'
  # wildcard1 / wildcard2 have a 8 chars match between read 0 and read 1;
  # the current version of readjoiner eliminates the wildcard-containing
  # sequences; update the test if this is changed
end

Name "gt readjoiner overlap: self-match"
Keywords "gt_readjoiner gt_readjoiner_overlap"
Test do
  # direct
  run_prefilter("#{$testdata}/readjoiner/self_spm.fas")
  run_overlap(8)
  grep last_stdout, 'number of irreducible suffix-prefix matches = 1'
  # with reverse complement of self
  run_prefilter("#{$testdata}/readjoiner/with_rc.fas")
  run_overlap(8)
  grep last_stdout, 'number of irreducible suffix-prefix matches = 1'
end

# the next test compares direct construction with complete graph + trans.red.
Name "gt readjoiner: transitive spm determination test - 1"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/transred_1.fas", "")
  run_overlap(4, "-singlestrand")
  run_assembly("-lengthcutoff 1 -depthcutoff 1")
  run_test "grep -f #{$testdata}/readjoiner/transred_1_targetseq.fas "+
           " reads.contigs.fas"
  run_overlap(4, "-singlestrand -elimtrans false")
  run_assembly("-redtrans -lengthcutoff 1 -depthcutoff 1")
  run_test "grep -f #{$testdata}/readjoiner/transred_1_targetseq.fas "+
           " reads.contigs.fas"
end

Name "gt readjoiner: transitive spm determination test - 2"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/transred_2.fas", "")
  run_overlap(4)
  run_assembly("-lengthcutoff 1 -depthcutoff 1")
  run_test "grep -f #{$testdata}/readjoiner/transred_1_targetseq.fas "+
           " reads.contigs.fas"
  run_overlap(4, "-elimtrans false")
  run_assembly("-redtrans -lengthcutoff 1 -depthcutoff 1")
  run_test "grep -f #{$testdata}/readjoiner/transred_1_targetseq.fas "+
           " reads.contigs.fas"
end

# currently the spm list obtained from transred_3 appears to be incorrect
# (the readset is varlen, currently not completely supported)
# the expected spm list from -l 10 -singlestrand is in transred_3.l10.ss.spm

Name "gt readjoiner: transitive spm determination test - 3"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/test_1.fas")
  run_overlap(39)
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(39, "-elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end


Name "gt readjoiner: transitive spm determination test - 4"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/test_2.fas")
  run_overlap(20)
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(20, "-elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end

Name "gt readjoiner: transitive spm determination test - 5"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/test_3.fas")
  run_overlap(20)
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(20, "-elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end

Name "gt readjoiner: transitive spm determination test - 6"
Keywords "gt_readjoiner"
Test do
  # this describes the next test case:
  #
  # 0 ctgataagtcccaggacttcagaagagctgtgag...ggacttcagaagagctgtgag
  # 2              ggacttcagaagagctgtgag...ggacttcagaagagctgtgagggtatggggacgg
  # 1          cccaggacttcagaagagctgtgag...ggacttcagaagagctgtgagggtatgggg
  #   |---A---|-B-|----------C----------|D|---------C----------|---E----|--F-|
  # 2                                      ggacttcagaagagctgtgagaccttggccaagtc..
  #                                       |---------C----------|-------D------..
  # assembly:
  # read_0:  ABCDC
  # read_1:   BCDCE
  # read_2:    CDCEF
  #
  # alternative assembly:
  # read_0:  ABCDC
  # read_1:   BCDCE
  # +
  # read_0: ABCDC
  # read_2:     CDCEF
  #
  #     ->2
  #    // ^
  #   0   |   only the longest of the two overlaps 0..2 is really transitive
  #    \  |
  #     ->1
  #

  run_prefilter("#{$testdata}/readjoiner/trans_and_submax_ovl.fas")
  run_overlap(20, "-singlestrand")
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(20, "-singlestrand -elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end

Name "gt readjoiner: transitive spm determination test - 7"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/large_count.fas")
  run_overlap(4)
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(4, "-elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end

Name "gt readjoiner: transitive spm determination test - 8"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/large_wset.fas")
  run_overlap(4)
  run_assembly
  run "mv reads.contigs.fas contigs"
  run_overlap(4, "-elimtrans false")
  run_assembly("-redtrans")
  run "diff reads.contigs.fas contigs"
end

Name "gt readjoiner: test different read lengths"
Keywords "gt_readjoiner"
Test do
  run_prefilter("#{$testdata}/readjoiner/70x_161nt.fas")
  run_overlap(30)
  run_assembly
end

# spmtest
[true, false].each do |mirrored|
  %w{contained_eqlen contained_varlen 30x_800nt 70x_100nt}.each do |fasta|
    Name "gt readjoiner spmtest#{' singlestrand' if !mirrored}: #{fasta}"
    Keywords "gt readjoiner gt_readjoiner_spmtest"
    Test do
      spmtest += " -singlestrand" if !mirrored
      encode_reads("#$testdata/readjoiner/#{fasta}.fas")
      run_test "#{spmtest} -test bruteforce -readset reads -l 32", :retval => 0
      bf = last_stdout
      run_test "#{spmtest} -test kmp -readset reads -l 32", :retval => 0
      kmp = last_stdout
      run "diff #{kmp} #{bf}"
      prepare_esa("reads", mirrored)
      if (!%w{70x_100nt}.include?(fasta))
        run_test "#{spmtest} -test gusfield -readset reads -l 32", :retval => 0
        run "sort -u #{last_stdout}"
        gf = last_stdout
        run "sort -u #{kmp}"
        kmp = last_stdout
        run "diff #{gf} #{kmp}"
      end
      if (!%w{contained_eqlen contained_varlen 30x_800nt}.include?(fasta))
        rdjO += " -singlestrand" if !mirrored
        run_test "#{rdjO} -readset reads -l 32 -showspm -elimtrans false -q",
          :retval => 0
        run "sort -u #{last_stdout}"
        our = last_stdout
        run "sort -u #{kmp}"
        kmp = last_stdout
        run "diff #{our} #{kmp}"
      end
    end
  end
end

# cnttest

[true, false].each do |mirrored|
  %w{contained_eqlen contained_varlen}.each do |fasta|
    Name "gt readjoiner cnttest#{' singlestrand' if !mirrored}: #{fasta}"
    Keywords "gt_readjoiner gt_readjoiner_cnttest"
    Test do
      cnttest += " -singlestrand" if !mirrored
      encode_reads("#$testdata/readjoiner/#{fasta}.fas")
      run_test "#{cnttest} -test bruteforce -readset reads", :retval => 0
      bf = last_stdout
      run_test "#{cnttest} -test kmp -readset reads", :retval => 0
      kmp = last_stdout
      run "diff #{kmp} #{bf}"
      prepare_esa("reads", mirrored)
      run_test "#{cnttest} -test esa -readset reads", :retval => 0
      esa = last_stdout
      run "diff #{esa} #{$testdata}readjoiner/#{fasta}#{'_ss' if !mirrored}.cnt"
      run "sort -u #{esa}"
      esas = last_stdout
      run "sort -u #{kmp}"
      kmps = last_stdout
      run "diff #{esas} #{kmps}"
    end
  end
end

def assert_nofspm(irreducible, transitive = nil)
  grep last_stdout,
    "number of irreducible suffix-prefix matches = #{irreducible}"
  if !transitive.nil?
    grep last_stdout,
      "number of transitive suffix-prefix matches = #{transitive}"
  end
end

[1, 2].each do |varlen_test|
  Name "gt encseq2spm: different minlen values - test #{varlen_test}"
  Keywords "gt_readjoiner gt_encseq2spm"
  Test do
    encode_reads("#$testdata/readjoiner/varlen_#{varlen_test}.fas")
    [10, 13, 15, 33, 65, 84].each do |minlen|
      run "#$bin/gt encseq2spm -l #{minlen} -ii reads"
    end
  end
end

Name "gt encseq2spm: different minlen values - test 3"
Keywords "gt_readjoiner gt_encseq2spm"
Test do
  run_prefilter("#$testdata/readjoiner/minlen_test.fas")
  run "#$bin/gt encseq2spm -l 14 -ii reads"
  run "#$bin/gt encseq2spm -l 15 -ii reads"
end

Name "gt readjoiner overlap: different min match lengths"
Keywords "gt_readjoiner gt_readjoiner_overlap"
Test do
  run_prefilter("#$testdata/readjoiner/minlen_test.fas")
  2.upto(40) do |minlen|
    run_overlap(minlen)
    assert_nofspm(40 - minlen)
  end
end

def unpack_uint_array_file(i_filename, o_filename, low_filter = 0)
  is64bit = Kernel.system("#{$bin}gt -64bit")
  sizeofint = is64bit ? 8 : 4
  unpackstr = is64bit ? "q" : "l"
  i_file = File.open(i_filename)
  o_file = File.open(o_filename, "w")
  loop do
    v = i_file.read(sizeofint)
    break if v.nil?
    value = v.unpack(unpackstr)[0]
    if value >= low_filter
      o_file.puts value
    end
  end
  i_file.close
  o_file.close
end

[true, false].each do |singlestrand|
  [0, 1000].each do |offset|
    %w{70x_100nt 30x_800nt}.each do |dataset|
      Name "gt readjoiner radixsort_str test "+
        "(#{dataset}, #{singlestrand ? 'single strand' : 'mirrored'}"+
        ", offset=#{offset})"
      Keywords "gt_readjoiner radixsort_str"
      Test do
        db = "#$testdata/readjoiner/#{dataset}.fas"
        run_prefilter(db, "-testrs -encseq no -q"+
                      "#{' -singlestrand' if singlestrand}"+
                      " -testrs-offset #{offset}")
        run "tail +#{offset+1} #{last_stdout}"
        radixsort_results = last_stdout
        run "#{$bin}gt suffixerator -suf -db #{db}"+
          "#{' -mirrored' unless singlestrand} -indexname i"
        unpack_uint_array_file("i.suf", "i.suf.txt", offset)
        run "diff i.suf.txt #{radixsort_results}"
      end
    end
  end
end

if $gttestdata

  # compare results with precalculated known results
  [700, 7000, 70000].each do |nofreads|
    readset="#{$gttestdata}/readjoiner/#{nofreads}x_100nt_reads"
    if File.exists?(readset)
      Name "gt readjoiner: #{nofreads}x100"
      Keywords "gt_readjoiner"
      Test do
        readset="#{$gttestdata}/readjoiner/#{nofreads}x_100nt_reads"
        run_prefilter(readset, "-cnt")
        run "#{cnttest} -test showlist -readset reads"
        run "diff #{last_stdout} #{readset}.cnt"
        run_overlap(45)
        run "#{spmtest} -test showlist -readset reads.0"
        run "sort #{last_stdout}"
        spm = last_stdout
        run "sort #{readset}.l45.spm"
        ref = last_stdout
        run "diff #{spm} #{ref}"
        run_assembly
        run "diff reads.contigs.fas #{readset}.l45.contigs"
      end
    end

    [161, 200, 300, 400, 600, 800, 1000].each do |len|
      reads = "#{$gttestdata}/readjoiner/#{nofreads}x_#{len}nt_reads"
      if File.exists?(reads)
        Name "gt readjoiner: #{nofreads}x#{len}"
        Keywords "gt_readjoiner"
        Test do
          run_prefilter(reads)
          run_overlap(30)
          run_assembly
        end
      end
    end
  end

  Name "gt readjoiner: multiple input files (eqlen)"
  Keywords "gt_readjoiner"
  Test do
    readset = "#{$gttestdata}/readjoiner/7000x_100nt_reads"
    run "head -n 3000 #{readset}"
    reads1 = last_stdout
    run "head -n 7000 #{readset} | tail -n 4000"
    reads2 = last_stdout
    run "head -n 13000 #{readset} | tail -n 6000"
    reads3 = last_stdout
    run "tail -n 1000 #{readset}"
    reads4 = last_stdout
    run_prefilter("#{reads1} #{reads2} #{reads3} #{reads4}")
    run_overlap(45)
    run "#{spmtest} -test showlist -readset reads.0"
    run "sort #{last_stdout}"
    spm = last_stdout
    run "sort #{readset}.l45.spm"
    ref = last_stdout
    run "diff #{spm} #{ref}"
  end

end

