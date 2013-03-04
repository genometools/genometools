Name "gt featureindex (empty file)"
Keywords "gt_featureindex"
Test do
  run "#{$bin}gt mkfeatureindex -filename tmp.db #{$testdata}/gt_view_prob_1.gff3"
  run "#{$bin}gt featureindex -filename tmp.db", :retval => 1
  grep(last_stderr, /no sequence regions in index/)
end

Name "gt featureindex (empty region)"
Keywords "gt_featureindex"
Test do
  run "#{$bin}gt mkfeatureindex -filename tmp.db #{$testdata}/gt_view_prob_2.gff3"
  run "#{$bin}gt featureindex -filename tmp.db"
  run "diff #{last_stdout} #{$testdata}/gt_view_prob_2.gff3"
end

Name "gt featureindex (parse error in GFF3)"
Keywords "gt_featureindex"
Test do
  run "#{$bin}gt mkfeatureindex -filename tmp.db #{$testdata}/gt_gff3_fail_1.gff3", :retval => 1
  grep(last_stderr, /has already been defined/)
  run "#{$bin}gt featureindex -filename tmp.db", :retval => 1
  grep(last_stderr, /no sequence regions in index/)
end

Name "gt featureindex (invalid sequence ID)"
Keywords "gt_featureindex"
Test do
  run "#{$bin}gt mkfeatureindex -filename tmp.db #{$testdata}/standard_gene_simple.gff3"
  run "#{$bin}gt featureindex -seqid foo -filename tmp.db", :retval => 1
  grep(last_stderr, /not exist/)
end

Name "gt featureindex (corrupt file)"
Keywords "gt_featureindex"
Test do
  File.open("corrupt.db", "w") do |file|
    file.write("sdfnhsnl")
  end
  run "#{$bin}gt featureindex -filename corrupt.db", :retval => 1
end

FEATUREINDEX_TEST_FILES = ["#{$testdata}/eden.gff3",
                           "#{$testdata}/standard_gene_simple.gff3",
                           "#{$testdata}/standard_gene_as_tree.gff3",
                           "#{$testdata}/standard_gene_with_introns_as_tree.gff3",
                           "#{$testdata}/encode_known_genes_Mar07.gff3"
                           ]

FEATUREINDEX_TEST_FILES.each do |file|
  Name "gt featureindex db vs. parser (#{File.basename(file)})"
  Keywords "gt_featureindex"
  Test do
    run "#{$bin}gt seqids #{file}"
    seqids = File.open(last_stdout).readlines
    run "#{$bin}gt mkfeatureindex -filename tmp.db #{file}", :maxtime => 1200
    seqids.each do |seqid|
      seqid.chomp!
      run "#{$bin}gt featureindex -seqid #{seqid} -retain no -filename tmp.db > out.gff3"
      run "#{$bin}gt gff3 -retainids no #{file} | #{$bin}gt select -seqid #{seqid}"
      run "diff out.gff3 #{last_stdout}"
    end
  end
end
