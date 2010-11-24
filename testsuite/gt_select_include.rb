Name "gt select test (no filter)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-seqid ctg123)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -seqid ctg123 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-seqid undef)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -seqid undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}header.gff3"
end

Name "gt select test (-source .)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -source . #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-source undef)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -source undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.header"
end

Name "gt select test (-maxgenelength)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenelength 8001 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenelength)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenelength 8000 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-mingenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -mingenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-mingenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -mingenescore .6 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-maxgenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenescore .4 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 1 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end
Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -strand + #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -strand - #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test("#{$bin}gt select -strand foo #{$testdata}standard_gene_as_tree.gff3",
           :retval => 1)
  grep last_stderr, /must be one of/
end

Name "gt select test (-overlap)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -overlap 2000 3000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-overlap)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -overlap 9001 10000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-minaveragessp)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -minaveragessp 0.5 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{last_stdout} #{$testdata}splice_site_prob.out"
end

Name "gt select test (-minaveragessp)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -minaveragessp 0.35 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{last_stdout} #{$testdata}splice_site_prob.gff3"
end

Name "gt select test (-hascds)"
Keywords "gt_select"
Test do
  run_test("#{$bin}gt select -hascds " +
           "#{$testdata}encode_known_genes_Mar07.gff3 | " +
           "#{$memcheck} #{$bin}gt stat", :maxtime => 120)
  run "diff #{last_stdout} #{$testdata}gt_select_encode.out"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1000 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1000-9000"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1001 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1001-9000"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1000 8999 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1000-8999"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1500000 1600000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}header.gff3"
end

Name "gt select test (-targetstrand)"
Keywords "gt_select targetstrand"
Test do
  run_test "#{$bin}gt select -targetstrand - #{$testdata}U89959_sas.gff3"
  run "diff #{last_stdout} #{$testdata}U89959_sas.minus_targets"
end

Name "gt select test (-targetstrand)"
Keywords "gt_select targetstrand"
Test do
  run_test "#{$bin}gt select -targetstrand + " +
            "#{$testdata}target_attribute_without_strand.gff3"
  run "diff #{last_stdout} #{$testdata}target_attribute_without_strand.gff3"
end

Name "gt select test (-targetbest, simple)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_simple_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_targetbest_simple_test.out"
end

Name "gt select test (-targetbest, complex)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_complex_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_targetbest_complex_test.out"
end

Name "gt select test (-targetbest, corrupt file)"
Keywords "gt_select targetbest"
Test do
  run_test("#{$bin}gt select -targetbest #{$testdata}corrupt_large.gff3",
           :retval => 1)
  grep last_stderr, "not a valid character"
end

Name "gt select test (-targetbest, multiple targets)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
  run      "diff #{last_stdout} " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
end

Name "gt select test (-rule_files, one file, node type in gff3)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_nodetype.lua -- " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-rule_files, one file, node type not in gff3)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_wrong_nodetype.lua -- " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-rule_files, two files, logic = AND)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_nodetype.lua " +
           "#{$testdata}gtscripts/filter_test_wrong_nodetype.lua -- " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-rule_files, two files, logic = OR)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_logic OR -rule_files " + 
           "#{$testdata}gtscripts/filter_test_nodetype.lua " +
           "#{$testdata}gtscripts/filter_test_wrong_nodetype.lua -- " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-rule_files, one file, wrong function name)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_wrong_function_name.lua -- " +
           "#{$testdata}standard_gene_as_tree.gff3", :retval => 1
  grep last_stderr, /error/
end

Name "gt select test (reading_frame_length % 3 != 0)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_orflength.lua -- " +
           "#{$testdata}filter_luafilter_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_luafilter_filtered_orfs.gff3"
end
  
Name "gt select test (check for LTR_retrotransposon and LTRs)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_LTR.lua -- " +
           "#{$testdata}filter_luafilter_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_luafilter_filtered_LTR.gff3"
end

Name "gt select test (min two orfs on forward strand)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_orf_pos_strand.lua -- " +
           "#{$testdata}filter_luafilter_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_luafilter_filtered_orf_pos.gff3"
end

Name "gt select test (orfs without frame attribute)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_frame_attribute.lua -- " +
           "#{$testdata}filter_luafilter_test_no_frame_attribute.gff3"
  run "diff #{last_stdout} #{$testdata}filter_luafilter_filtered_orf_frame.gff3"
end

Name "gt select test (dropped to file (1))"
 Keywords "gt_select"
 Test do
   run_test "#{$bin}gt select -dropped_file nh_file01.gff3 -rule_files " + 
            "#{$testdata}gtscripts/filter_test_orflength.lua -- " +
            "#{$testdata}filter_luafilter_test.gff3"
   run "diff nh_file01.gff3 #{$testdata}filter_nh_file01.gff3"
end
  
Name "gt select test (dropped to file (2))"
 Keywords "gt_select"
 Test do
   run_test "#{$bin}gt select -dropped_file nh_file02.gff3 -rule_files " + 
            "#{$testdata}gtscripts/filter_test_LTR.lua -- " +
            "#{$testdata}filter_luafilter_test.gff3"
   run "diff nh_file02.gff3 #{$testdata}filter_nh_file02.gff3"
end

Name "gt select test (dropped to file (3))"
 Keywords "gt_select"
 Test do
   run_test "#{$bin}gt select -dropped_file nh_file03.gff3 -rule_files " + 
            "#{$testdata}gtscripts/filter_test_orf_pos_strand.lua -- " +
            "#{$testdata}filter_luafilter_test.gff3"
   run "diff nh_file03.gff3 #{$testdata}filter_nh_file03.gff3"
end

Name "gt select test (dropped to file (4))"
 Keywords "gt_select"
 Test do
   run_test "#{$bin}gt select -dropped_file nh_file04.gff3 -rule_files " + 
            "#{$testdata}gtscripts/filter_test_frame_attribute.lua -- " +
            "#{$testdata}filter_luafilter_test_no_frame_attribute.gff3"
   run "diff nh_file04.gff3 #{$testdata}filter_nh_file04.gff3"
end

Name "gt select test (lua syntax fail)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -rule_files " + 
           "#{$testdata}gtscripts/filter_test_syntax_fail.lua -- " +
           "#{$testdata}filter_luafilter_test_no_frame_attribute.gff3",
           :retval => 1
  grep last_stderr, /error/
end
