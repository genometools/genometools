Name "gt sketch short test"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png #{$testdata}gff3_file_1_short.txt", \
           :maxtime => 600
  run "test -e out.png"
end

Name "gt sketch short test (stdin)"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png < #{$testdata}gff3_file_1_short.txt", \
           :maxtime => 600
  run "test -e out.png"
end

Name "gt sketch short test (unknown output format)"
Keywords "gt_sketch"
Test do
  run_test("#{$bin}gt sketch -format unknown out.png " + \
           "#{$testdata}gff3_file_1_short.txt", :retval => 1, :maxtime => 600)
  grep(last_stderr, /must be one of:/)
end

Name "gt sketch short test (unwriteable PNG file)"
Keywords "gt_sketch"
Test do
  run "touch unwriteable.png"
  run "chmod u-w unwriteable.png"
  run_test("#{$bin}gt sketch -force unwriteable.png " + \
           "#{$testdata}gff3_file_1_short.txt", :retval => 1, :maxtime => 600)
  grep(last_stderr, /an I\/O error occurred/)
end

Name "gt sketch short test (unwriteable PDF file)"
Keywords "gt_sketch"
Test do
  run "touch unwriteable.pdf"
  run "chmod u-w unwriteable.pdf"
  run_test("#{$bin}gt sketch -format pdf -force unwriteable.pdf " + \
           "#{$testdata}gff3_file_1_short.txt", :retval => 1, :maxtime => 600)
  grep(last_stderr, /Permission denied/)
end

Name "gt sketch short test (nonexistant style file)"
Keywords "gt_sketch"
Test do
  run_test("#{$bin}gt sketch -style foo.bar out.png " + \
           "#{$testdata}gff3_file_1_short.txt", :retval => 1, :maxtime => 600)
  grep(last_stderr, /style file 'foo.bar' does not exist/)
end

Name "gt sketch short test (invalid style file)"
Keywords "gt_sketch"
Test do
  run "echo thisisnotlua > foo.bar"
  run_test("#{$bin}gt sketch  -style foo.bar out.png " + \
           "#{$testdata}gff3_file_1_short.txt", :retval => 1, :maxtime => 600)
  grep(last_stderr, /cannot run style file/)
end

Name "gt sketch prob 1"
Keywords "gt_sketch"
Test do
  run_test("#{$bin}gt sketch out.png #{$testdata}gt_view_prob_1.gff3", \
           :retval => 1, :maxtime => 600)
end

Name "gt sketch prob 2"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png #{$testdata}gt_view_prob_2.gff3", \
           :maxtime => 600
  run "test -e out.png"
end

Name "gt sketch track IDs to style file"
Keywords "gt_sketch"
Test do
  run "cp #{$testdata}eden.gff3 eden2.gff3"
  run_test "#{$bin}gt sketch -force -style " + \
           "#{$testdata}trackname1.style out.png #{$testdata}eden.gff3 " + \
           "eden2.gff3"
  run "diff #{last_stdout} #{$testdata}trackname1.out"
  run "test -e out.png"
end

Name "gt sketch multiline without parent"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png #{$testdata}gt_sketch_multiline_without_parent.gff3", :maxtime => 600
  run "test -e out.png"
end

Name "gt sketch pipe"
Keywords "gt_sketch"
Test do
  run "#{$bin}gt gff3 #{$testdata}gff3_file_1_short.txt > in.gff3"
  run_test "#{$bin}gt sketch -pipe out.png in.gff3 > out.gff3", :maxtime => 600
  run "diff in.gff3 out.gff3"
end

Name "gt sketch streams <-> file output"
Keywords "gt_sketch streams annotationsketch"
Test do
  m = `#{$bin}gt sketch -help`.match(/graphics format\s+choose from ([^ ]+)/m)
  raise TestFailedError if m.nil?
  m[1].chomp.split("|").each do |format|
    run_test "#{$bin}gt sketch -format #{format} out.#{format} #{$testdata}eden.gff3"
    run_test "#{$bin}gt sketch -streams -format #{format} streamout.#{format} #{$testdata}eden.gff3"
    # some formats will diff in their creation date, remove it
    run "sed -i '/CreationDate/d' out.#{format} streamout.#{format}"
    run "diff out.#{format} streamout.#{format}"
  end
end

Name "gt sketch -showrecmaps"
Keywords "gt_sketch showrecmaps"
Test do
  run_test "#{$bin}gt sketch -showrecmaps out.png " +
           "#{$testdata}standard_gene_as_tree.gff3", :maxtime => 600
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.recmaps"
end

Name "gt sketch -showrecmaps (normal text size)"
Keywords "gt_sketch showrecmaps"
Test do
  run_test "#{$bin}gt sketch -showrecmaps out.png " +
           "#{$testdata}gt_sketch_textwidth.gff3", :maxtime => 600
  run "diff #{last_stdout} #{$testdata}gt_sketch_textwidth_0.recmaps"
end

Name "gt sketch -showrecmaps (narrow image)"
Keywords "gt_sketch showrecmaps"
Test do
  run_test "#{$bin}gt sketch -width 300 -showrecmaps out.png " +
           "#{$testdata}gt_sketch_textwidth.gff3", :maxtime => 600
  run "diff #{last_stdout} #{$testdata}gt_sketch_textwidth_1.recmaps"
end

Name "gt sketch -showrecmaps (large text size)"
Keywords "gt_sketch showrecmaps"
Test do
  run_test "#{$bin}gt sketch -style #{$testdata}bigfonts.style " + \
           "-showrecmaps out.png #{$testdata}gt_sketch_textwidth.gff3", \
           :maxtime => 600
  run "diff #{last_stdout} #{$testdata}gt_sketch_textwidth_2.recmaps"
end

Name "gt sketch for transcript (neg. coords in ruler)"
Keywords "gt_sketch annotationsketch neg_coords"
Test do
  run_test "#{$bin}gt sketch " + \
           "-style #{$testdata}transcript.style " + \
           "transcript.png " + \
           "#{$testdata}transcript.gff3 ", \
           :maxtime => 600
end

Name "gt sketch runtime Lua failures"
Keywords "gt_sketch lua"
Test do
  Dir.glob("#{$testdata}fail*style") do |file|
    run_test("#{$bin}gt sketch -style #{file} " + \
             "out.png #{$testdata}eden.gff3", \
             :maxtime => 600, :retval => 1)
    grep(last_stderr, 'attempt to call global \'fail\'')
  end
end

Name "sketch_constructed (C)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}examples/sketch_constructed " + \
           "#{$cur}/gtdata/sketch/default.style sketch_constructed.png", \
           :maxtime => 600
end

Name "sketch_parsed (C)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}examples/sketch_parsed " + \
           "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
           "#{$testdata}standard_gene_with_introns_as_tree.gff3 ", \
           :maxtime => 600
end

Name "sketch_parsed reverse order (C)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}examples/sketch_parsed_with_ordering " + \
           "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
           "#{$testdata}eden.gff3 ", \
           :maxtime => 600
  run "diff #{last_stdout} #{$testdata}order_sketch_out.txt"
end

Name "sketch_constructed (Lua)"
Keywords "gt_sketch gt_scripts annotationsketch"
Test do
  run_test "#{$bin}gt #{$cur}/gtscripts/sketch_constructed.lua " + \
           "#{$cur}/gtdata/sketch/default.style sketch_constructed.png", \
           :maxtime => 600
end

Name "sketch_parsed (Lua)"
Keywords "gt_sketch gt_scripts annotationsketch"
Test do
  run_test "#{$bin}gt #{$cur}/gtscripts/sketch_parsed.lua " + \
           "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
           "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
           :maxtime => 600
end

if python_tests_runnable? and not $arguments["nocairo"] then
  Name "sketch_constructed (Python)"
  Keywords "gt_sketch gt_python annotationsketch"
  Test do
    run_python "#{$cur}/gtpython/sketch_constructed.py " + \
               "#{$cur}/gtdata/sketch/default.style sketch_constructed.png", \
               :maxtime => 600
  end

  Name "sketch_parsed (Python)"
  Keywords "gt_sketch gt_python annotationsketch"
  Test do
    run_python "#{$cur}/gtpython/sketch_parsed.py " + \
               "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
               "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
               :maxtime => 600
  end

  Name "sketch_parsed reverse order (Python)"
  Keywords "gt_sketch gt_python annotationsketch"
  Test do
    run_python "#{$testdata}gtpython/sketch_parsed_with_ordering.py " + \
               "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
               "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
               :maxtime => 600
  end

  Name "sketch_parsed invalid order (Python, string != int)"
  Keywords "gt_sketch gt_python annotationsketch"
  Test do
    run_python "#{$testdata}gtpython/sketch_parsed_with_invalid_ordering.py " + \
               "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
               "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
               :maxtime => 600
    grep(last_stderr, /Track ordering function must return a number/)
  end

  Name "sketch_parsed invalid order (Python, None)"
  Keywords "gt_sketch gt_python annotationsketch"
  Test do
    run_python "#{$testdata}gtpython/sketch_parsed_with_invalid_ordering_2.py " + \
               "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
               "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
               :maxtime => 600
    grep(last_stderr, /Track ordering function must return a number/)
  end

  Name "Python runtime style failures"
    Keywords "gt_sketch gt_python annotationsketch"
    Test do
    Dir.glob("#{$testdata}fail*style") do |file|
      run_python "#{$cur}/gtpython/sketch_parsed.py " + \
                 "#{file} sketch_parsed.png " + \
                 "#{$testdata}eden.gff3", \
                 :maxtime => 600, :retval => 1
      grep(last_stderr, 'attempt to call global \'fail\'')

      run_python "#{$cur}/gtpython/sketch_constructed.py " + \
                 "#{file} sketch_constructed.png", \
                 :maxtime => 600, :retval => 1
      grep(last_stderr, 'attempt to call global \'fail\'')

      run_python "#{$testdata}gtpython/style_serialize.py #{file}", :retval => 1
      grep(last_stderr, 'expected boolean, number, or string')
    end
  end
end

if ruby_tests_runnable? and not $arguments["nocairo"] then
  Name "sketch_constructed (Ruby)"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    run_ruby "#{$cur}/gtruby/sketch_constructed.rb " + \
             "#{$cur}/gtdata/sketch/default.style sketch_constructed.png", \
             :maxtime => 600
  end

  Name "sketch_parsed (Ruby)"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    run_ruby "#{$cur}/gtruby/sketch_parsed.rb " + \
             "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
             "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
             :maxtime => 600
  end

  Name "sketch_parsed reverse order (Ruby)"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    run_ruby "#{$testdata}gtruby/sketch_parsed_with_ordering.rb " + \
             "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
             "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
             :maxtime => 600
  end

  Name "sketch_parsed invalid order (Ruby, string != numeric)"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    run_ruby "#{$testdata}gtruby/sketch_parsed_with_invalid_ordering.rb " + \
             "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
             "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
             :maxtime => 600, :retval => 1
    grep(last_stderr, /Track ordering callback must return a number/)
  end

  Name "sketch_parsed invalid order (Ruby, nil)"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    run_ruby "#{$testdata}gtruby/sketch_parsed_with_invalid_ordering_2.rb " + \
             "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " + \
             "#{$testdata}standard_gene_with_introns_as_tree.gff3", \
             :maxtime => 600, :retval => 1
    grep(last_stderr, /Track ordering callback must return a number/)
  end

  Name "Ruby runtime style failures"
  Keywords "gt_sketch gt_ruby annotationsketch"
  Test do
    Dir.glob("#{$testdata}fail*style") do |file|
      run_ruby "#{$cur}/gtruby/sketch_parsed.rb " + \
               "#{file} sketch_parsed.png " + \
               "#{$testdata}eden.gff3", \
               :maxtime => 600, :retval => 1
      grep(last_stderr, 'attempt to call global \'fail\'')

      run_ruby "#{$cur}/gtruby/sketch_constructed.rb " + \
               "#{file} sketch_constructed.png", \
               :maxtime => 600, :retval => 1
      grep(last_stderr, 'attempt to call global \'fail\'')

      run_ruby "#{$testdata}gtruby/style_serialize.rb #{file}", :retval => 1
      grep(last_stderr, 'expected boolean, number, or string')
    end
  end
end
