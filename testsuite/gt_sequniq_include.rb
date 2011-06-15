["", " -r", " -seqit", " -seqit -r"].each do |opt|
  Name "gt sequniq#{opt} 2xfoo test"
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} #{$testdata}foofoo.fas"
    run "diff #{$last_stdout} #{$testdata}foo.fas"
  end

  Name "gt sequniq#{opt} 3xfoo test "
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} #{$testdata}foofoofoo.fas"
    run "diff #{$last_stdout} #{$testdata}foo.fas"
  end
end

["", " -seqit"].each do |opt|
  Name "gt sequniq#{opt} foo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} #{$testdata}foorcfoo.fas"
    run "diff #{$last_stdout} #{$testdata}foorcfoo.fas"
  end

  Name "gt sequniq#{opt} -r foo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} -r #{$testdata}foorcfoo.fas"
    run "diff #{$last_stdout} #{$testdata}foo.fas"
  end

  Name "gt sequniq#{opt} 2xfoo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} #{$testdata}foorcfoofoo.fas"
    run "diff #{$last_stdout} #{$testdata}foorcfoo.fas"
  end

  Name "gt sequniq#{opt} -r 2xfoo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    run_test "#{$bin}gt sequniq#{opt} -r #{$testdata}foorcfoofoo.fas"
    run "diff #{$last_stdout} #{$testdata}foo.fas"
  end
end