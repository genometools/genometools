Name "gt exercise affinealign test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise affinealign #{$testdata}gt_affinealign_test1_seq1.fas #{$testdata}gt_affinealign_test1_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_affinealign_test1.out"
end

Name "gt exercise affinealign test 2 (default)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise affinealign #{$testdata}gt_affinealign_test2_seq1.fas #{$testdata}gt_affinealign_test2_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_affinealign_test2.out1"
end

Name "gt exercise affinealign test 2 (no gap opening costs)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise affinealign -gapopen 0 #{$testdata}gt_affinealign_test2_seq1.fas #{$testdata}gt_affinealign_test2_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_affinealign_test2.out2"
end

Name "gt exercise align test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise align #{$testdata}gt_align_test.in1 #{$testdata}gt_align_test.in2"
  run "diff #{$last_stdout} #{$testdata}gt_align_test.out"
end

Name "gt exercise align -all test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise align -all #{$testdata}gt_align_test.in1 #{$testdata}gt_align_test.in2"
  run "diff #{$last_stdout} #{$testdata}gt_align_test.all"
end

Name "gt exercise align (qgramdist seqs)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise align #{$testdata}gt_qgramdist_test_seq1.fas #{$testdata}gt_qgramdist_test_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_qgramdist_test_align.out"
end

Name "gt exercise assemblegreedy 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy " +
           "#{$testdata}gt_assemblegreedy_fragments_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_1.assembly"
end

Name "gt exercise assemblegreedy 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy " +
           "#{$testdata}gt_assemblegreedy_fragments_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_2.assembly"
end

Name "gt exercise assemblegreedy 3"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy " +
           "#{$testdata}gt_assemblegreedy_fragments_3.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_3.assembly"
end

Name "gt exercise assemblegreedy 2 (-minlength)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy -minlength 6 " +
           "#{$testdata}gt_assemblegreedy_fragments_2.fas"
  run "diff #{$last_stdout} " +
       "#{$testdata}gt_assemblegreedy_fragments_2.assembly_minlength6"
end

Name "gt exercise assemblegreedy 1 (-showoverlaps)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy -showoverlaps " +
           "#{$testdata}gt_assemblegreedy_fragments_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_1.overlaps"
end

Name "gt exercise assemblegreedy 2 (-showoverlaps)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy -showoverlaps " +
           "#{$testdata}gt_assemblegreedy_fragments_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_2.overlaps"
end

Name "gt exercise assemblegreedy 1 (-showpath)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy -showpath " +
           "#{$testdata}gt_assemblegreedy_fragments_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_1.path"
end

Name "gt exercise assemblegreedy 2 (-showpath)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise assemblegreedy -showpath " +
           "#{$testdata}gt_assemblegreedy_fragments_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_assemblegreedy_fragments_2.path"
end

Name "gt exercise assemblegreedy nonexistent file"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise assemblegreedy #{$testdata}nonexistent_file",
           :retval => 1)
end

# XXX: fix this test
=begin
Name "gt exercise assemblegreedy large test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt shredder -overlap 20 #{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt exercise assemblegreedy #{$last_stdout}"
  run_test "#{$bin}gt bioseq -showfasta -width 70 #{$last_stdout}"
  run "tail -n 1529 #{$last_stdout} > assembly.txt"
  run "tail -n 1529 #{$testdata}U89959_genomic.fas > original.txt"
  run "diff -u original.txt assembly.txt"
end
=end

Name "gt exercise blastenv (script)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise blastenv #{$testdata}identity_score_matrix " +
           "GGCCGC"
  run "egrep '^[GC]{4}' #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gt_blastenv_test_script.out"
end

Name "gt exercise blastenv (lookahead)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise blastenv -k 18 #{$testdata}BLOSUM62 " +
           "SIEQSIEQSIEQSIEQSIEQSIEQSIEQSIEQSIEQSIEQSQSQSQSSSQQSQSS" +
           "QSQSSSQSSQSQSSQSQQSQSQSSQS"
  run "diff #{$last_stdout} #{$testdata}gt_blastenv_test_lookahead.out"
end

Name "gt exercise casino test "
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise casino `cat #{$testdata}t_casino.in`"
  run "diff #{$last_stdout} #{$testdata}t_casino.out"
end

Name "gt exercise coin test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise coin `cat #{$testdata}t_coin_1.in`"
  run "diff #{$last_stdout} #{$testdata}t_coin_1.out"
end

Name "gt exercise coin test 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise coin `cat #{$testdata}t_coin_2.in`"
  run "diff #{$last_stdout} #{$testdata}t_coin_2.out"
end

1.upto(3) do |i|
  Name "gt exercise fastaparser success #{i}"
  Keywords "gt_exercise"
  Test do
    run_test "#{$bin}gt exercise fastaparser " +
             "#{$testdata}gt_bioseq_succ_#{i}.fas"
  end
end

1.upto(7) do |i|
  Name "gt exercise fastaparser failure #{i}"
  Keywords "gt_exercise"
  Test do
    run_test("#{$bin}gt exercise fastaparser " +
             "#{$testdata}gt_bioseq_fail_#{i}.fas", :retval => 1)
  end
end

if $gttestdata then # XXX: hack to avoid testing matchcount during normal runs
                    # (is quite expensive due to amount of calls)
  Name "gt exercise matchcount test"
  Keywords "gt_exercise"
  Test do
    1.upto(5) do |k|
      for u in ['a', 'b', 'ab', 'bb', 'aba', 'aab', 'baa', 'abba', 'baaba', \
                'aaaba', 'bbbbb']
        for v in ['a', 'b', 'aa', 'ab', 'ba', 'aba', 'abb', 'aaaa', 'baaba', \
                  'bbbbb', 'aaaba']
          run_test "#{$bin}gt exercise matchcount #{k} #{u} #{v}"
          run "cat #{$last_stdout} >> result.txt"
        end
      end
    end
    run "cmp result.txt #{$testdata}gt_matchcount_test.out"
  end
end

Name "gt exercise neighborjoining test example"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise neighborjoining example"
  # XXX: the output of the following command is not the same on all platforms
  # (probably due to floating point differences)
  # run "diff #{$last_stdout} #{$testdata}gt_neighborjoining_example.out"
end

Name "gt exercise nussinov_rna_fold test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise nussinov_rna_fold 3 -3 -2 -1 ACCCU"
  run "diff #{$last_stdout} #{$testdata}gt_nussinov_rna_fold_test_1.out"
end

Name "gt exercise nussinov_rna_fold test 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise nussinov_rna_fold 3 -3 -2 -1 GGGAAAUCC"
  run "diff #{$last_stdout} #{$testdata}gt_nussinov_rna_fold_test_2.out"
end

Name "gt exercise qgramdist test (q=2)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise qgramdist -q 2 #{$testdata}gt_qgramdist_test_seq1.fas #{$testdata}gt_qgramdist_test_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_qgramdist_test_q2.out"
end

Name "gt exercise qgramdist test (q=3)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise qgramdist -q 3 #{$testdata}gt_qgramdist_test_seq1.fas #{$testdata}gt_qgramdist_test_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_qgramdist_test_q3.out"
end

Name "gt exercise linearalign test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise linearalign #{$testdata}gt_align_test.in1 #{$testdata}gt_align_test.in2"
  run "diff #{$last_stdout} #{$testdata}gt_linearalign_test.out"
end

Name "gt exercise nussinov_rna_fold test (tRNA glutamine)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise nussinov_rna_fold 3 -3 -2 -1 ggttccatggtgtaatggttagcactctggactctgaatccagcgatccgagttcaaatctcggtggaacct"
  run "diff #{$last_stdout} #{$testdata}trna_glutamine.out"
end

Name "gt exercise nussinov_rna_fold test (script example)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise nussinov_rna_fold 1 -1 -1 -1 GGGAAAUCC"
  run "diff #{$last_stdout} #{$testdata}script_example.out"
end

Name "gt exercise markovchain"
Keywords "gt_exercise gt_markovchain"
Test do
  run_test "#{$bin}gt exercise markovchain " +
           "#{$testdata}gt_markovchain_example.txt bACGTe"
  run "diff #{$last_stdout} #{$testdata}gt_markovchain_example.out"
end

Name "gt exercise msaparse fail 1"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise msaparse #{$testdata}gt_msaparse_fail_1.fas",
           :retval => 1)
end

Name "gt exercise msaparse fail 2"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise msaparse #{$testdata}gt_msaparse_fail_2.fas",
           :retval => 1)
end

Name "gt exercise msaparse test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msaparse -show #{$testdata}gt_msaparse_test_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_msaparse_test_1.fas"
end

Name "gt exercise msaparse test 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msaparse -show #{$testdata}gt_msaparse_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_msaparse_test_2.fas"
end

Name "gt exercise msaparse test 2 (consensus)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msaparse -consensus #{$testdata}gt_msaparse_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_msaparse_test_2.cons"
end

Name "gt exercise msaparse test 2 (sum of pairs)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msaparse -sumofpairs #{$testdata}gt_msaparse_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_msaparse_test_2.sop"
end

Name "gt exercise msmatch test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msmatch `cat #{$testdata}test_1_multiset.txt` `cat #{$testdata}test_1_text.txt`"
  run "diff #{$last_stdout} #{$testdata}test_1_result.txt"
end

Name "gt exercise msmatch test 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msmatch `cat #{$testdata}test_2_multiset.txt` `cat #{$testdata}test_2_text.txt`"
  run "diff #{$last_stdout} #{$testdata}test_2_result.txt"
end

Name "gt exercise msmatch test 3"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise msmatch `cat #{$testdata}test_3_multiset.txt` `cat #{$testdata}test_3_text.txt`"
  run "diff #{$last_stdout} #{$testdata}test_3_result.txt"
end

Name "gt exercise multilcp -help"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise multilcp -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt exercise multilcp -noop"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise multilcp -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

Name "gt exercise multilcp test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise multilcp aacacac ccacacaa"
  run "diff #{$last_stdout} #{$testdata}gt_multilcp_test.out"
end

Name "gt exercise multilcp fail"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise multilcp '' ''", :retval => 1)
  grep $last_stderr, "sequence of length 0 not allowed"
end

Name "gt exercise scorefasta test"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise scorefasta -q 1 AGTGCACACATC ATCACACTTAGC"
  run "diff #{$last_stdout} #{$testdata}score_6.txt"

  run_test "#{$bin}gt exercise scorefasta -q 2 AGTGCACACATC ATCACACTTAGC"
  run "diff #{$last_stdout} #{$testdata}score_4.txt"

  run_test "#{$bin}gt exercise scorefasta -q 3 AGTGCACACATC ATCACACTTAGC"
  run "diff #{$last_stdout} #{$testdata}score_3.txt"

  run_test "#{$bin}gt exercise scorefasta -q 2 AGCGATAG AGTGACAG"
  run "diff #{$last_stdout} #{$testdata}score_3.txt"

  run_test "#{$bin}gt exercise scorefasta -q 3 AGCGATAG AGTGACAG"
  run "diff #{$last_stdout} #{$testdata}score_0.txt"

  run_test "#{$bin}gt exercise scorefasta -q 3 AG AGTGACAG"
  run "diff #{$last_stdout} #{$testdata}score_0.txt"

  run_test "#{$bin}gt exercise scorefasta -q 3 AGCGATAG AG"
  run "diff #{$last_stdout} #{$testdata}score_0.txt"
end

Name "gt exercise scorematrix test BLOSUM62"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise scorematrix #{$testdata}BLOSUM62"
  run "diff #{$last_stdout} #{$testdata}BLOSUM62.out"
end

Name "gt exercise scorematrix test simple scorematrix"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise scorematrix #{$testdata}simple_scorematrix.in"
  run "diff #{$last_stdout} #{$testdata}simple_scorematrix.out"
end

Name "gt exercise swalign test (empty)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign #{$testdata}BLOSUM62 #{$testdata}gt_swalign_empty_seq1.fas #{$testdata}gt_swalign_empty_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt exercise swalign test (simple)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign #{$testdata}BLOSUM62 #{$testdata}gt_swalign_simple_seq1.fas #{$testdata}gt_swalign_simple_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_simple.out"
end

Name "gt exercise swalign test (simple self 1)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign #{$testdata}BLOSUM62 #{$testdata}gt_swalign_simple_seq1.fas #{$testdata}gt_swalign_simple_seq1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_simple.out_s1"
end

Name "gt exercise swalign test (simple self 2)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign #{$testdata}BLOSUM62 #{$testdata}gt_swalign_simple_seq2.fas #{$testdata}gt_swalign_simple_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_simple.out_s2"
end

Name "gt exercise swalign test (script)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign -indelscore -1 #{$testdata}simple_scorematrix.in #{$testdata}gt_swalign_script_seq1.fas #{$testdata}gt_swalign_script_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_script.out"
end

Name "gt exercise swalign test (script self 1)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign -indelscore -1 #{$testdata}simple_scorematrix.in #{$testdata}gt_swalign_script_seq1.fas #{$testdata}gt_swalign_script_seq1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_script.out_s1"
end

Name "gt exercise swalign test (script self 2)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise swalign -indelscore -1 #{$testdata}simple_scorematrix.in #{$testdata}gt_swalign_script_seq2.fas #{$testdata}gt_swalign_script_seq2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_swalign_script.out_s2"
end

Name "gt exercise scorematrix test unsymmetric scorematrix"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_1",
           :retval => 1)
end

Name "gt exercise scorematrix test double char on alphabet line"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_2",
           :retval => 1)
end

Name "gt exercise scorematrix test double char on score line"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_3",
           :retval => 1)
end

Name "gt exercise scorematrix test illegal token on alpha line 1"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_4",
           :retval => 1)
end

Name "gt exercise scorematrix test illegal token on alpha line 2"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_5",
           :retval => 1)
end

Name "gt exercise scorematrix test illegal token on score line"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}corrupt_scorematrix_6",
           :retval => 1)
end

Name "gt exercise scorematrix test empty file as scorematrix"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}empty_file",
           :retval => 1)
end

Name "gt exercise scorematrix test directory as scorematrix"
Keywords "gt_exercise"
Test do
  run_test("#{$bin}gt exercise scorematrix #{$testdata}",
           :retval => 1)
end

Name "gt exercise translate"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise translate `cat #{$testdata}gt_translate.in`"
  run "diff #{$last_stdout} #{$testdata}gt_translate.out"
end

Name "gt exercise translate (too short)"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise translate gg"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt exercise upgma test example"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise upgma example"
  run "diff #{$last_stdout} #{$testdata}gt_upgma_example.out"
end

Name "gt exercise upgma test 1"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise upgma #{$testdata}gt_upgma_test_1.fas"
  run "diff #{$last_stdout} #{$testdata}gt_upgma_test_1.out"
end

Name "gt exercise upgma test 2"
Keywords "gt_exercise"
Test do
  run_test "#{$bin}gt exercise upgma #{$testdata}gt_upgma_test_2.fas"
  run "diff #{$last_stdout} #{$testdata}gt_upgma_test_2.out"
end

i = 1
for infile in `ls #{$testdata}/consensus_sa/*.in` do
  Name "gt exercise consensus_sa test #{i}"
  Keywords "gt_exercise gt_consensus_sa"
  Test do
    run_test "#{$bin}gt exercise consensus_sa #{infile}"
    outfile = infile.gsub(/\.in$/, ".out")
    run "diff #{$last_stdout} #{outfile}"
  end
  i += 1
end
