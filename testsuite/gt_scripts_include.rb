Name "arg passing test 1"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/arg.lua"
  grep($last_stdout, /^arg\[0\]=.*gt.*arg\.lua/)
end

Name "arg passing test 2"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/arg.lua foo"
  grep($last_stdout, /^arg\[0\]=.*gt.*arg\.lua/)
  grep($last_stdout, /^arg\[1\]=foo$/)
end

Name "bittab bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/bittab.lua"
end

1.upto(4) do |i|
  Name "csa_stream bindings #{i}"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt  #{$testdata}/gtscripts/csa_stream.lua " +
             "#{$testdata}/gt_csa_prob_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}/gt_csa_prob_#{i}.out"
  end
end

1.upto(14) do |i|
  Name "cds_stream bindings #{i}"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt  #{$testdata}/gtscripts/cds_stream.lua " +
             "#{$testdata}/gt_cds_test_#{i}.fas " +
             "#{$testdata}/gt_cds_test_#{i}.in"
    run "sed 's/gtscript/gt cds/' #{$last_stdout}"
    run "diff #{$last_stdout} #{$testdata}/gt_cds_test_#{i}.out"
  end
end

Name "genome_node bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_node.lua"
end

Name "genome_stream bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_stream.lua #{$testdata}"
end

Name "genome_stream bindings (output stream)"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/gff3.lua #{$testdata}/gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "genome_visitor bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_visitor.lua #{$testdata}/gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "range bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/range.lua"
end

Name "memleak.lua"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}gtscripts/memleak.lua #{$testdata}standard_gene_as_tree.gff3"
end

Name "scorematrix2c"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/scorematrix2c.lua #{$testdata}/BLOSUM62"
  run "diff #{$last_stdout} #{$testdata}blosum62.c"
end

Name "require 'gtlua'"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/require_gtlua.lua"
end

Name "extract_swalign"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/extract_swalign.lua #{$cur} " +
           "#{$testdata}BLOSUM62"
  run "cd swalign && ${MAKE:-make}"
  if not File.exists?("swalign.tar.gz") then
    raise TestFailed, "file 'swalign.tar.gz' does not exist"
  end
end

Name "extract_swalign -sol"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/extract_swalign.lua -sol " +
           "#{$cur} #{$testdata}BLOSUM62"
  run "cd swalign && ${MAKE:-make} test"
  if not File.exists?("swalign.tar.gz") then
    raise TestFailed, "file 'swalign.tar.gz' does not exist"
  end
end

Name "extract_affinealign"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/extract_affinealign.lua #{$cur}"
  run "cd affinealign && ${MAKE:-make}"
  if not File.exists?("affinealign.tar.gz") then
    raise TestFailed, "file 'affinealign.tar.gz' does not exist"
  end
end

Name "extract_linearalign"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/extract_linearalign.lua #{$cur}"
  run "cd linearalign && ${MAKE:-make} test"
  if not File.exists?("linearalign.tar.gz") then
    raise TestFailed, "file 'linearalign.tar.gz' does not exist"
  end
end

Name "extract_assemblegreedy"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/extract_assemblegreedy.lua " +
           "#{$cur}"
  run "cd assemblegreedy && ${MAKE:-make} test"
  if not File.exists?("assemblegreedy.tar.gz") then
    raise TestFailed, "file 'assemblegreedy.tar.gz' does not exist"
  end
end

Name "extract_assemblegreedy_scaffold"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt " +
           "#{$testdata}../gtscripts/extract_assemblegreedy_scaffold.lua " +
           "#{$cur}"
  run "cd assemblegreedy && ${MAKE:-make}"
  if not File.exists?("assemblegreedy.tar.gz") then
    raise TestFailed, "file 'assemblegreedy.tar.gz' does not exist"
  end
end

Name "extract_consensus_sa"
Keywords "gt_scripts extract"
Test do
  run_test "#{$bin}gt " +
           "#{$testdata}../gtscripts/extract_consensus_sa.lua #{$cur}"
  run "cd consensus_sa && ${MAKE:-make}"
  if not File.exists?("consensus_sa.tar.gz") then
    raise TestFailed, "file 'consensus_sa.tar.gz' does not exist"
  end
end

Name "LPeg library"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$cur}/src/external/lpeg-0.7/test.lua"
end

Name "MD5 library"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$cur}/src/external/md5-1.1.2/tests/test.lua"
end

Name "gtdoc"
Keywords "gt_scripts gtdoc"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/gtdoc.lua #{$cur}"
end

Name "gtdoc -v"
Keywords "gt_scripts gtdoc"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/gtdoc.lua -v #{$cur}"
end

Name "gtdoc -html"
Keywords "gt_scripts gtdoc"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/gtdoc.lua -html #{$cur}"
end

Name "gtdoc -html -v"
Keywords "gt_scripts gtdoc"
Test do
  run_test "#{$bin}gt #{$testdata}../gtscripts/gtdoc.lua -html -v #{$cur}"
end

if $arguments["libgtview"] then
  Name "feature_index and feature_stream bindings"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}/gtscripts/feature_stuff.lua #{$testdata}"
    run "env LC_ALL=C sort #{$last_stdout}"
    run "grep -v '^##sequence-region' #{$testdata}gff3_file_1_short_sorted.txt | diff #{$last_stdout} -"
  end

  Name "libgtview bindings"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}/gtscripts/view.lua test.png #{$testdata}gff3_file_1_short.txt"
  end

  Name "show_seqids"
  Keywords "gt_scripts"
  Test do
    run_test("#{$bin}gt #{$testdata}/gtscripts/show_seqids.lua #{$testdata}encode_known_genes_Mar07.gff3", :maxtime => 100)
    run "diff #{$last_stdout} #{$testdata}encode_known_genes_Mar07.seqids"
  end

  Name "evalviz.lua test 1"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}../gtscripts/evalviz.lua png_files #{$testdata}/gt_eval_test_1.in #{$testdata}/gt_eval_test_1.in"
    run "grep -v seqid #{$last_stdout}"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_1.out"
  end

=begin XXX: takes too long
  Name "evalviz.lua test 2"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}../gtscripts/evalviz.lua png_files #{$testdata}/gt_evalviz_test.reality #{$testdata}/gt_evalviz_test.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_evalviz_test.out"
  end
=end
end
