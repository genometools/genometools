["", " -seqit", " -seqit -r"].each do |opt|
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

Name "gt sequniq -seqit foo + rc(foo) test "
Keywords "gt_sequniq"
Test do
  run_test "#{$bin}gt sequniq -seqit #{$testdata}foorcfoo.fas"
  run "diff #{$last_stdout} #{$testdata}foorcfoo.fas"
end

Name "gt sequniq -seqit -r foo + rc(foo) test "
Keywords "gt_sequniq"
Test do
  run_test "#{$bin}gt sequniq -seqit -r #{$testdata}foorcfoo.fas"
  run "diff #{$last_stdout} #{$testdata}foo.fas"
end

Name "gt sequniq -seqit 2xfoo + rc(foo) test "
Keywords "gt_sequniq"
Test do
  run_test "#{$bin}gt sequniq -seqit #{$testdata}foorcfoofoo.fas"
  run "diff #{$last_stdout} #{$testdata}foorcfoo.fas"
end

Name "gt sequniq -seqit -r 2xfoo + rc(foo) test "
Keywords "gt_sequniq"
Test do
  run_test "#{$bin}gt sequniq -seqit -r #{$testdata}foorcfoofoo.fas"
  run "diff #{$last_stdout} #{$testdata}foo.fas"
end
