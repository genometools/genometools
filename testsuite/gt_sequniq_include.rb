require "fileutils"

["", " -rev", " -seqit", " -seqit -rev"].each do |opt|
  Name "gt sequniq#{opt} 2xfoo test"
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foofoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} foofoo.fas"
    run "diff #{last_stdout} #{$testdata}foo.fas"
  end

  Name "gt sequniq#{opt} 3xfoo test "
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foofoofoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} foofoofoo.fas"
    run "diff #{last_stdout} #{$testdata}foo.fas"
  end
end

["", " -seqit"].each do |opt|
  Name "gt sequniq#{opt} foo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foorcfoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} foorcfoo.fas"
    run "diff #{last_stdout} #{$testdata}foorcfoo.fas"
  end

  Name "gt sequniq#{opt} -rev foo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foorcfoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} -rev foorcfoo.fas"
    run "diff #{last_stdout} #{$testdata}foo.fas"
  end

  Name "gt sequniq#{opt} 2xfoo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foorcfoofoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} foorcfoofoo.fas"
    run "diff #{last_stdout} #{$testdata}foorcfoo.fas"
  end

  Name "gt sequniq#{opt} -rev 2xfoo + rc(foo) test "
  Keywords "gt_sequniq"
  Test do
    FileUtils.copy("#{$testdata}foorcfoofoo.fas", ".")
    run_test "#{$bin}gt sequniq#{opt} -rev foorcfoofoo.fas"
    run "diff #{last_stdout} #{$testdata}foo.fas"
  end
end

Name "gt sequniq (revbug, no -rev)"
Keywords "gt_sequniq"
Test do
  FileUtils.copy("#{$testdata}gt_sequniq_rev_bug.fas", ".")
  run_test "#{$bin}gt sequniq gt_sequniq_rev_bug.fas"
  run "diff #{last_stdout} #{$testdata}gt_sequniq_rev_bug.fas"
end

Name "gt sequniq (revbug, -rev)"
Keywords "gt_sequniq"
Test do
  FileUtils.copy("#{$testdata}gt_sequniq_rev_bug.fas", ".")
  run_test "#{$bin}gt sequniq -rev gt_sequniq_rev_bug.fas"
  run "diff #{last_stdout} #{$testdata}gt_sequniq_rev_bug.out"
end
