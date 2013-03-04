Name "gtpython: genome_stream bindings (output stream)"
Keywords "gt_python"
Test do
  run_python "#{$testdata}gtpython/gff3.py #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gtpython: genome_visitor bindings (output stream)"
Keywords "gt_python"
Test do
  run_python "#{$testdata}gtpython/genome_visitor.py #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gtpython: feature_index and feature_stream bindings"
Keywords "gt_python"
Test do
  run_python "#{$testdata}gtpython/feature_stuff.py " +
             "#{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{last_stdout}"
  run "grep -v '^##sequence-region' #{$testdata}gff3_file_1_short_sorted.txt | diff #{last_stdout} -"
end

if not $arguments["nocairo"] then
  Name "gtpython: AnnotationSketch bindings (valid gff3 file)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/sketch.py test.png " +
               "#{$testdata}gff3_file_1_short.txt"
  end

  Name "gtpython: AnnotationSketch bindings (corrupt gff3 file)"
  Keywords "gt_python"
  Test do
    run_python("#{$testdata}gtpython/sketch.py test.png #{$testdata}corrupt.gff3",
               :retval => 1)
    grep last_stderr, "GenomeTools error"
  end

  Name "gtpython: AnnotationSketch bindings (nonexistent gff3 file)"
  Keywords "gt_python"
  Test do
    run_python("#{$testdata}gtpython/sketch.py test.png " +
               "#{$testdata}nonexistent_file", :retval => 1)
    grep last_stderr, "GenomeTools error"
  end

  Name "gtpython: AnnotationSketch bindings (simple sketch)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/sketch_simple.py test.png " +
               "#{$testdata}gff3_file_1_short.txt"
  end

  Name "gtpython: AnnotationSketch bindings (PNG stream)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/sketch_stream.py test.png " +
               "#{$testdata}gff3_file_1_short.txt"
  end

  Name "gtpython: AnnotationSketch bindings (TrackSelectorFunc)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/block_stuff.py " +
             "#{$testdata}gff3_file_1_short.txt"
    run "env LC_ALL=C sort #{last_stdout}"
    run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.blocks"
  end

  Name "gtpython: AnnotationSketch bindings (style)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/style.py #{$cur}/gtdata/sketch/default.style"
  end

  Name "gtpython: AnnotationSketch bindings (error reporting)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/sketch-failures.py " +
               "#{$testdata}gff3_file_1_short.txt"
  end

  Name "gtpython: AnnotationSketch bindings (Graphics)"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/graphics_stuff.py " +
             "#{$testdata}graphics_curve_test_coords.txt " +
             "out.svg"
    # will fail e.g. if cairo toy font setup is different from test machine
    # disabled for now
    # run "diff out.svg #{$testdata}graphics_test.out"
  end

  Name "gtpython: AnnotationSketch bindings (FeatureNode(Iterator))"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/feature_node.py"
  end

  Name "gtpython: show_seqids"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/show_seqids.py #{$testdata}encode_known_genes_Mar07.gff3"
    run "diff #{last_stdout} #{$testdata}encode_known_genes_Mar07.seqids"
  end

  Name "gtpython: used_types"
  Keywords "gt_python"
  Test do
    run_python "#{$testdata}gtpython/used_types.py " +
               "#{$testdata}standard_gene_as_tree.gff3"
    run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.types"
  end

  Name "gtpython: show_recmaps"
  Keywords "gt_python showrecmaps"
  Test do
    run_python "#{$testdata}gtpython/show_recmaps.py " +
               "#{$testdata}standard_gene_as_tree.gff3"
    run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.hotspots"
  end

  Name "gtpython: unicode strings"
  Keywords "gt_python"
  Test do
    run_python("#{$testdata}gtpython/unicode_strings.py #{$cur}/gtdata/sketch/default.style test.png")
  end
end


Name "gtpython: unittests"
Keywords "gt_python unittests"
Test do
  run_python "#{$gtpython}/tests/__init__.py "
  grep last_stderr, "OK"
end
