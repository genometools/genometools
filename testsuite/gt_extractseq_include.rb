Name "gt extractseq -help"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt extractseq -noop"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

Name "gt extractseq test stdin"
Keywords "gt_extractseq"
Test do
  run "cat #{$testdata}foo.fas | #{$memcheck} #{$bin}gt extractseq -match foo"
  run "diff #{$last_stdout} #{$testdata}foo.fas"
end

Name "gt extractseq test foo"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match foo #{$testdata}foo.fas"
  run "diff #{$last_stdout} #{$testdata}foo.fas"
end

Name "gt extractseq test foo width 4"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match foo -width 4 #{$testdata}foo.fas"
  run "diff #{$last_stdout} #{$testdata}foo_width4.fas"
end

Name "gt extractseq test bar"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match bar -width 4 #{$testdata}bar.fas"
  run "diff #{$last_stdout} #{$testdata}bar.fas"
end

Name "gt extractseq test baz"
Keywords "gt_extractseq"
Test do
  run "cat #{$testdata}foo.fas | #{$memcheck} #{$bin}gt extractseq -match baz - #{$testdata}bar.fas"
  grep($last_stdout, ".", true)
end

Name "gt extractseq test foo|bar"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match 'foo|bar' #{$testdata}foo.fas #{$testdata}bar.fas"
  run "diff #{$last_stdout} #{$testdata}foobar.fas"
end

Name "gt extractseq test '(foo'"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -match '(foo' #{$testdata}foo.fas", :retval => 1)
end

Name "gt extractseq test corrupt"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -match foo #{$testdata}corrupt.fas", :retval => 1)
  grep $last_stderr, "first character"
end

Name "gt extractseq test corrupt (stdin)"
Keywords "gt_extractseq"
Test do
  run("cat #{$testdata}corrupt.fas | #{$memcheck} #{$bin}gt extractseq -match foo",
      :retval => 1)
  grep $last_stderr, "first character"
end
