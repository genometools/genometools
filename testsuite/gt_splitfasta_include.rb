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

Name "gt splitfasta (-splitdesc, nonexistent file)"
Keywords "gt_splitfasta"
Test do
  run_test("#{$bin}gt splitfasta -splitdesc . #{$testdata}nonexistent_file",
           :retval => 1)
  grep $last_stderr, /does not exist or is not readable/
end

Name "gt splitfasta (empty file)"
Keywords "gt_splitfasta"
Test do
  run_test("#{$bin}gt splitfasta #{$testdata}empty_file", :retval => 1)
  grep $last_stderr, /is empty/
end


Name "gt splitfasta (corrupt file)"
Keywords "gt_splitfasta"
Test do
  run_test("#{$bin}gt splitfasta #{$testdata}gt_bioseq_fail_2.fas",
           :retval => 1)
  grep $last_stderr, /file is not in FASTA format/
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

if $gttestdata then
  Name "gt splitfasta (large data, compressed)"
  Keywords "gt_splitfasta"
  Test do
    run "cp #{$gttestdata}ltrharvest/s_cer/chrAll_before-1997-10-01.fsa.gz " +
        "test.fas.gz"
    run_test("#{$bin}gt splitfasta -targetsize 2 test.fas.gz", :maxtime => 160)
    # check sample
    if not File.exists?("test.fas.2.gz") then
      raise TestFailed, "file 'test.fas.2.gz' does not exist"
    end
  end

  Name "gt splitfasta (large data, uncompressed)"
  Keywords "gt_splitfasta"
  Test do
    run "cp #{$gttestdata}ltrharvest/s_cer/chrAll_before-1997-10-01.fsa.gz " +
        "test.fas.gz"
    run "gunzip test.fas.gz"
    run_test "#{$bin}gt splitfasta -targetsize 2 test.fas"
    # check sample
    if not File.exists?("test.fas.2") then
      raise TestFailed, "file 'test.fas.2' does not exist"
    end
  end

  Name "gt splitfasta (large data, file exists)"
  Keywords "gt_splitfasta"
  Test do
    run "cp #{$gttestdata}ltrharvest/s_cer/chrAll_before-1997-10-01.fsa.gz " +
        "test.fas.gz"
    run "touch test.fas.2.gz"
    run_test("#{$bin}gt splitfasta -targetsize 1 test.fas.gz",
             :retval => 1, :maxtime => 120)
    grep $last_stderr, /exists already/
  end
end
