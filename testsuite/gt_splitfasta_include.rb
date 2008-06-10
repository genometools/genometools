Name "gt splitfasta (default)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}U89959_genomic.fas ."
  run_test "#{$bin}gt splitfasta U89959_genomic.fas"
  if not File.exists?("U89959_genomic.fas.1") then
    raise TestFailed, "file 'U89959_genomic.fas.1' does not exist"
  end
end

Name "gt splitfasta (compressed file)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}U89959_genomic.fas ."
  run "gzip U89959_genomic.fas"
  run_test "#{$bin}gt splitfasta U89959_genomic.fas.gz"
  if not File.exists?("U89959_genomic.fas.1.gz") then
    raise TestFailed, "file 'U89959_genomic.fas.1.gz' does not exist"
  end
end

Name "gt splitfasta (file exists)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}U89959_genomic.fas ."
  run "touch U89959_genomic.fas.1"
  run_test("#{$bin}gt splitfasta U89959_genomic.fas", :retval => 1)
  grep $last_stderr, /exists already/
end

Name "gt splitfasta (-force)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}U89959_genomic.fas ."
  run "touch U89959_genomic.fas.1"
  run_test "#{$bin}gt splitfasta -force U89959_genomic.fas"
  if not File.exists?("U89959_genomic.fas.1") then
    raise TestFailed, "file 'U89959_genomic.fas.1' does not exist"
  end
end

Name "gt splitfasta (-targetsize)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}U89959_genomic.fas ."
  run_test "#{$bin}gt splitfasta -targetsize 1 U89959_genomic.fas"
  if not File.exists?("U89959_genomic.fas.1") then
    raise TestFailed, "file 'U89959_genomic.fas.1' does not exist"
  end
end

Name "gt splitfasta (-splitdesc)"
Keywords "gt_splitfasta"
Test do
  run_test "#{$bin}gt splitfasta -splitdesc . #{$testdata}foobar.fas"
  if not File.exists?("foo.fas") then
    raise TestFailed, "file 'foo.fas' does not exist"
  end
  if not File.exists?("bar.fas") then
    raise TestFailed, "file 'bar.fas' does not exist"
  end
end

Name "gt splitfasta (-splitdesc, file exists)"
Keywords "gt_splitfasta"
Test do
  run "touch foo.fas"
  run_test("#{$bin}gt splitfasta -splitdesc . #{$testdata}foobar.fas",
           :retval => 1)
  grep $last_stderr, /exists already/
end

Name "gt splitfasta (-splitdesc -force)"
Keywords "gt_splitfasta"
Test do
  run "touch foo.fas"
  run_test "#{$bin}gt splitfasta -force -splitdesc . #{$testdata}foobar.fas"
  if not File.exists?("foo.fas") then
    raise TestFailed, "file 'foo.fas' does not exist"
  end
  if not File.exists?("bar.fas") then
    raise TestFailed, "file 'bar.fas' does not exist"
  end
end
Name "gt splitfasta (-splitdesc, compressed file)"
Keywords "gt_splitfasta"
Test do
  run "cp #{$testdata}foobar.fas ."
  run "gzip foobar.fas"
  run_test "#{$bin}gt splitfasta -splitdesc . foobar.fas.gz"
  if not File.exists?("foo.fas.gz") then
    raise TestFailed, "file 'foo.fas.gz' does not exist"
  end
  if not File.exists?("bar.fas.gz") then
    raise TestFailed, "file 'bar.fas.gz' does not exist"
  end
end


