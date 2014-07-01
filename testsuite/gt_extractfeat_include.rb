Name "gt extractfeat -seqfile test 1"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}gt_extractfeat_succ_1.fas -matchdesc #{$testdata}gt_extractfeat_succ_1.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -seqfile test 1 (compressed)"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}gt_extractfeat_succ_1.fas.gz -matchdesc  #{$testdata}gt_extractfeat_succ_1.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -seqfile test 1 (with -o)"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -o test.fas -type gene -seqfile #{$testdata}gt_extractfeat_succ_1.fas.gz -matchdesc  #{$testdata}gt_extractfeat_succ_1.gff3"
  run "diff test.fas #{$testdata}gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -seqfile test 2"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}gt_extractfeat_succ_2.fas -matchdesc #{$testdata}gt_extractfeat_succ_2.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_2.out1"
end

Name "gt extractfeat -seqfile test 3"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -seqfile #{$testdata}gt_extractfeat_succ_2.fas -matchdesc #{$testdata}gt_extractfeat_succ_2.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_2.out2"
end

Name "gt extractfeat -seqfile test 4"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -join -seqfile #{$testdata}gt_extractfeat_succ_2.fas -matchdesc #{$testdata}gt_extractfeat_succ_2.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_2.out3"
end

Name "gt extractfeat -seqfile test 5"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -join -seqfile #{$testdata}gt_extractfeat_succ_3.fas -matchdesc #{$testdata}gt_extractfeat_succ_3.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_3.out"
end

Name "gt extractfeat -regionmapping fail 1 (no mapping file)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}nonexistent_file #{$testdata}gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep(last_stderr, "cannot run file");
end

Name "gt extractfeat -regionmapping fail 2 (empty file)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}empty_file #{$testdata}gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep(last_stderr, "'mapping' is not defined ");
end

Name "gt extractfeat -regionmapping fail 3 (wrong type)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}regionmapping_1.lua #{$testdata}gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep(last_stderr, "'mapping' must be either a table or a function ");
end

Name "gt extractfeat -regionmapping fail 3 (nil mapping)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}regionmapping_2.lua #{$testdata}gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep(last_stderr, "is nil");
end

Name "gt extractfeat -regionmapping fail 4 (non string mapping)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}regionmapping_3.lua #{$testdata}gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep(last_stderr, "is not a string");
end

Name "gt extractfeat -regionmapping fail 5 (function returns nil)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}regionmapping_5.lua #{$testdata}gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep(last_stderr, "function 'mapping' must return a string");
end

Name "gt extractfeat -regionmapping test 1 (mapping table)"
Keywords "gt_extractfeat"
Test do
  run "env GT_TESTDATA=#{$testdata} #{$memcheck} #{$bin}gt extractfeat -type gene -regionmapping #{$testdata}regionmapping_4.lua #{$testdata}gt_extractfeat_succ_1.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -regionmapping test 1 (mapping function)"
Keywords "gt_extractfeat"
Test do
  run "env GT_TESTDATA=#{$testdata} #{$memcheck} #{$bin}gt extractfeat -type gene -regionmapping #{$testdata}regionmapping_6.lua #{$testdata}gt_extractfeat_succ_1.gff3"
end

Name "gt extractfeat error message"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt gff3 -offset 1000 #{$testdata}gt_extractfeat_succ_1.gff3 | " +
      "#{$bin}gt extractfeat -type gene -seqfile " +
      "#{$testdata}gt_extractfeat_succ_1.fas -matchdesc -", :retval => 1
  grep last_stderr, "Has the sequence-region to sequence mapping been defined correctly"
end

Name "gt extractfeat -translate"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-matchdesc -type CDS -join -translate #{$testdata}U89959_cds.gff3"
  run "diff #{last_stdout} #{$testdata}U89959_cds.fas"
end

Name "gt extractfeat -translate (phases)"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}gt_extractfeat_phase.fas " +
      "-matchdesc -type CDS -join -translate #{$testdata}gt_extractfeat_phase_fix.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_phase_fix.out"
  run "#{$bin}gt extractfeat -seqfile #{$testdata}gt_extractfeat_phase.fas " +
      "-matchdesc -type CDS -join -translate #{$testdata}gt_extractfeat_phase.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_phase.out"
end

Name "gt extractfeat -help"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -help"
  grep(last_stdout, "Lua");
end

Name "gt extractfeat -seqid"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-type CDS -join -seqid " +
      "#{$testdata}gt_extractfeat_seqid_target.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_seqid.fas"
end

Name "gt extractfeat -target"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-type mRNA -join -target " +
      "#{$testdata}gt_extractfeat_seqid_target.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_target.fas"
end

Name "gt extractfeat -seqid -target"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-type CDS -join -seqid -target " +
      "#{$testdata}gt_extractfeat_seqid_target.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_seqid_target.fas"
end

Name "gt extractfeat -retainids"
Keywords "gt_extractfeat"
Test do
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-type CDS -retainids -join -translate " +
      "#{$testdata}gt_extractfeat_retainids.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_retainids_join.fas"
  run "#{$bin}gt extractfeat -seqfile #{$testdata}U89959_genomic.fas " +
      "-type CDS -retainids -translate " +
      "#{$testdata}gt_extractfeat_retainids.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_retainids.fas"
end

Name "gt extractfeat compare query method combinations"
Keywords "gt_extractfeat query_methods"
Test do
  run "#{$bin}gt encseq encode -lossless -indexname ./idx " +
      "#{$testdata}gt_extractfeat_mappings.fas"
  ["", ".md5"].each do |md5|
    ["-usedesc","-matchdesc"].each do |method|
      run "#{$bin}gt extractfeat #{method} " +
          "-seqfile #{$testdata}gt_extractfeat_mappings.fas " +
          "-type gene #{$testdata}gt_extractfeat_mappings#{md5}.gff3 "
      run "diff #{last_stdout} #{$testdata}gt_extractfeat_mappings_ref#{md5}.fas"
      run "#{$bin}gt extractfeat #{method} -seqfiles " +
          "#{$testdata}gt_extractfeat_mappings_sep1.fas " +
          "#{$testdata}gt_extractfeat_mappings_sep2.fas " +
          "-type gene #{$testdata}gt_extractfeat_mappings#{md5}.gff3 "
      run "diff #{last_stdout} #{$testdata}gt_extractfeat_mappings_ref#{md5}.fas"
      run "#{$bin}gt extractfeat #{method} -encseq idx " +
          "-type gene #{$testdata}gt_extractfeat_mappings#{md5}.gff3 "
      run "diff #{last_stdout} #{$testdata}gt_extractfeat_mappings_ref#{md5}.fas"
    end
    run "env GT_TESTDATA=#{$testdata} #{$memcheck} #{$bin}gt extractfeat " +
        "-regionmapping #{$testdata}gt_extractfeat_mappings_seprm.lua " +
        "-type gene #{$testdata}gt_extractfeat_mappings#{md5}.gff3 "
    run "diff #{last_stdout} #{$testdata}gt_extractfeat_mappings_ref#{md5}.fas"
  end
end

Name "gt extractfeat -matchdescstart"
Keywords "gt_extractfeat matchdescstart"
Test do
  ["-seqfiles","-seqfile"].each do |method|
    run "#{$bin}gt extractfeat #{method} #{$testdata}gt_extractfeat_matchdescstart_1.fas " +
        "-type gene -matchdesc #{$testdata}gt_extractfeat_matchdescstart_1.gff3", \
        :retval => 1
    grep(last_stderr, "could match more than one sequence")
    run "#{$bin}gt extractfeat #{method} #{$testdata}gt_extractfeat_matchdescstart_1.fas " +
        "-type gene -matchdescstart #{$testdata}gt_extractfeat_matchdescstart_1.gff3"
    run "diff #{last_stdout} #{$testdata}gt_extractfeat_matchdescstart_1.out"
    run "#{$bin}gt extractfeat #{method} #{$testdata}gt_extractfeat_matchdescstart_2.fas " +
        "-type gene -matchdescstart #{$testdata}gt_extractfeat_matchdescstart_1.gff3", \
        :retval => 1
    grep(last_stderr, "could match more than one sequence")
  end
  run "#{$bin}gt encseq encode -lossless -indexname foo " +
      "#{$testdata}gt_extractfeat_matchdescstart_1.fas"
  run "#{$bin}gt extractfeat -encseq foo " +
      "-type gene -matchdesc #{$testdata}gt_extractfeat_matchdescstart_1.gff3", \
      :retval => 1
  grep(last_stderr, "could match more than one sequence")
  run "#{$bin}gt extractfeat -encseq foo " +
      "-type gene -matchdescstart #{$testdata}gt_extractfeat_matchdescstart_1.gff3"
  run "diff #{last_stdout} #{$testdata}gt_extractfeat_matchdescstart_1.out"
  run "#{$bin}gt encseq encode -lossless -indexname foo " +
      "#{$testdata}gt_extractfeat_matchdescstart_2.fas"
  run "#{$bin}gt extractfeat -encseq foo " +
      "-type gene -matchdescstart #{$testdata}gt_extractfeat_matchdescstart_1.gff3", \
      :retval => 1
  grep(last_stderr, "could match more than one sequence")
end
