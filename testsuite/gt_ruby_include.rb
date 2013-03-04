rv = RUBY_VERSION.match(/^(\d+\.\d+)/)
if (!rv.nil? && rv[1].to_f < 1.9)  then
  Name "gtruby: genome_stream bindings (output stream)"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/gff3.rb #{$testdata}gff3_file_1_short.txt"
    run "env LC_ALL=C sort #{last_stdout}"
    run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
  end

  Name "gtruby: genome_visitor bindings (output stream)"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/genome_visitor.rb #{$testdata}gff3_file_1_short.txt"
    run "env LC_ALL=C sort #{last_stdout}"
    run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
  end

  Name "gtruby: feature_index and feature_stream bindings"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/feature_stuff.rb " +
             "#{$testdata}gff3_file_1_short.txt"
    run "env LC_ALL=C sort #{last_stdout}"
    run "grep -v '^##sequence-region' #{$testdata}gff3_file_1_short_sorted.txt | diff #{last_stdout} -"
  end
  if not $arguments["nocairo"] then
    Name "gtruby: AnnotationSketch bindings (valid gff3 file)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/sketch.rb test.png " +
               "#{$testdata}gff3_file_1_short.txt"
    end

    Name "gtruby: AnnotationSketch bindings (corrupt gff3 file)"
    Keywords "gt_ruby"
    Test do
      run_ruby("#{$testdata}gtruby/sketch.rb test.png #{$testdata}corrupt.gff3",
               :retval => 1)
      grep last_stderr, "GenomeTools error"
    end

    Name "gtruby: AnnotationSketch bindings (nonexistent gff3 file)"
    Keywords "gt_ruby"
    Test do
      run_ruby("#{$testdata}gtruby/sketch.rb test.png " +
               "#{$testdata}nonexistent_file", :retval => 1)
      grep last_stderr, "GenomeTools error"
    end

    Name "gtruby: AnnotationSketch bindings (PNG stream)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/sketch_stream.rb test.png " +
               "#{$testdata}gff3_file_1_short.txt"
    end

    Name "gtruby: AnnotationSketch bindings (TrackSelectorFunc)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/block_stuff.rb " +
               "#{$testdata}gff3_file_1_short.txt"
      run "env LC_ALL=C sort #{last_stdout}"
      run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.blocks"
    end

    Name "gtruby: AnnotationSketch bindings (style)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/style.rb #{$cur}/gtdata/sketch/default.style"
    end

    Name "gtruby: AnnotationSketch bindings (error reporting)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/sketch-failures.rb " +
               "#{$testdata}gff3_file_1_short.txt"
    end

    Name "gtruby: AnnotationSketch bindings (Graphics)"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/graphics_stuff.rb " +
               "#{$testdata}graphics_curve_test_coords.txt " +
               "out.svg"
      # will fail e.g. if cairo toy font setup is different from test machine
      # disabled for now
      # run "diff out.svg #{$testdata}graphics_test.out"
    end

    Name "gtruby: AnnotationSketch bindings (FeatureNode(Iterator))"
    Keywords "gt_ruby"
    Test do
      run_ruby "#{$testdata}gtruby/feature_node.rb"
    end

    Name "gtruby: show_recmaps"
    Keywords "gt_ruby showrecmaps"
    Test do
      run_ruby "#{$testdata}gtruby/show_recmaps.rb " +
               "#{$testdata}standard_gene_as_tree.gff3"
      run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.hotspots"
    end
  end

  Name "gtruby: Encseq bindings"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/encseq.rb"
  end

  Name "gtruby: show_seqids"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/show_seqids.rb #{$testdata}encode_known_genes_Mar07.gff3"
    run "diff #{last_stdout} #{$testdata}encode_known_genes_Mar07.seqids"
  end

  Name "gtruby: used_types"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/used_types.rb " +
             "#{$testdata}standard_gene_as_tree.gff3"
    run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.types"
  end

  Name "gtruby: {Comment,Sequence,Region,Meta,EOF}Node classes"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/node_types.rb"
  end

  Name "gtruby: CustomStream/CustomVisitor basic tests"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/custom_stuff.rb #{$testdata}eden.gff3"
    run "diff #{last_stdout} #{$testdata}custom_streams_ref.txt"
  end

  Name "gtruby: CustomStream/CustomVisitor all node types"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/custom_visitor.rb"
  end

  Name "gtruby: Range class"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/range.rb"
  end

  Name "gtruby: TypeChecker class"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/type_checker.rb #{$cur}/gtdata/obo_files/so.obo"
  end

  Name "gtruby: TypeChecker class (failure)"
  Keywords "gt_ruby"
  Test do
    run_ruby "#{$testdata}gtruby/type_checker.rb #{$gtdata}/obo_files/so", :retval => 1
  end
end
