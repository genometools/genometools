Name "gt extractseq -help"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -help"
  grep last_stdout, "Report bugs to"
end

Name "gt extractseq -noop"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -noop", :retval => 1)
  grep last_stderr, "unknown option"
end

Name "gt extractseq test stdin"
Keywords "gt_extractseq"
Test do
  run "cat #{$testdata}foo.fas | #{$memcheck} #{$bin}gt extractseq -match foo"
  run "diff #{last_stdout} #{$testdata}foo.fas"
end

Name "gt extractseq test foo"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match foo #{$testdata}foo.fas"
  run "diff #{last_stdout} #{$testdata}foo.fas"
end

Name "gt extractseq test foo width 4"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match foo -width 4 #{$testdata}foo.fas"
  run "diff #{last_stdout} #{$testdata}foo_width4.fas"
end

Name "gt extractseq test bar"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match bar -width 4 #{$testdata}bar.fas"
  run "diff #{last_stdout} #{$testdata}bar.fas"
end

Name "gt extractseq test baz"
Keywords "gt_extractseq"
Test do
  run "cat #{$testdata}foo.fas | #{$memcheck} #{$bin}gt extractseq -match baz - #{$testdata}bar.fas"
  grep(last_stdout, ".", true)
end

Name "gt extractseq test foo|bar"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -match 'foo|bar' #{$testdata}foo.fas #{$testdata}bar.fas"
  run "diff #{last_stdout} #{$testdata}foobar.fas"
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
  grep last_stderr, "cannot guess"
end

Name "gt extractseq test corrupt (stdin)"
Keywords "gt_extractseq"
Test do
  run("cat #{$testdata}corrupt.fas | #{$memcheck} #{$bin}gt extractseq -match foo",
      :retval => 1)
end

Name "gt extractseq -frompos"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -frompos 5 -topos 12 #{$testdata}foobar.fas"
  run "diff #{last_stdout} #{$testdata}frompos.fas"
end

Name "gt extractseq -frompos (stdin)"
Keywords "gt_extractseq"
Test do
  run "cat  #{$testdata}foobar.fas | #{$memcheck} #{$bin}gt extractseq " +
      "-frompos 5 -topos 12"
  run "diff #{last_stdout} #{$testdata}frompos.fas"
end

Name "gt extractseq -frompos (fail 1)"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -frompos 5 -topos 17 #{$testdata}foobar.fas",
           :retval => 1
  grep last_stderr, "larger than"
end

Name "gt extractseq -frompos (fail 2)"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -frompos 18 -topos 17 #{$testdata}foobar.fas",
           :retval => 1
  grep last_stderr, "must be <= argument"
end

Name "gt extractseq marker.fas"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -frompos 21 -topos 40 #{$testdata}marker.fas"
  run "diff #{last_stdout} #{$testdata}marker.out"
end

Name "gt extractseq -keys from fastafile U89959"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -keys #{$testdata}U89959_ginums.txt " +
           "#{$testdata}U89959_ests.fas"
  run "grep -v '^#' #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}U89959_ginums.out"
end

Name "gt extractseq -keys from fastafile at1MB"
Keywords "gt_extractseq"
Test do
  run "sed  -e '/^[^\\>]/d' -e 's/^>gi\|\\([^\|]*\\).*/\\1/' #{$testdata}at1MB"
  run_test "#{$bin}gt extractseq -keys #{last_stdout} -width 70 " +
           "#{$testdata}at1MB"
  run "grep -v '^#' #{last_stdout}"
  run "cmp #{last_stdout} #{$testdata}at1MB"
end

Name "gt extractseq -keys from fastafile TrEMBL"
Keywords "gt_extractseq"
Test do
  run_test "#{$bin}gt extractseq -keys #{$testdata}trembl-keys.txt -width 60 " +
           "#{$testdata}trembl.faa"
  run "grep -v '^#' #{last_stdout}"
  run "cmp #{last_stdout} #{$testdata}trembl.faa"
end

Name "gt extractseq -keys from fastafile (corrupt)"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -keys #{$testdata}U89959_ginums.corrupt " +
           "#{$testdata}U89959_ests.fas", :retval => 1)
end

Name "gt extractseq -keys from fastafile (fail)"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt extractseq -keys #{$testdata}U89959_ginums.txt",
           :retval => 1)
  grep last_stderr, /requires at least one file argument/
end

Name "gt suffixerator -kys with key of length 0"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt suffixerator -protein -ssp -tis -des -sds -kys " +
             "-db #{$testdata}trembl-emptykey.faa",:retval => 1)
end

Name "gt suffixerator -kys with keys of different length"
Keywords "gt_extractseq"
Test do
  run_test("#{$bin}gt suffixerator -protein -ssp -tis -des -sds -kys " +
             "-db #{$testdata}trembl-difflenkey.faa",:retval => 1)
end

if $gttestdata then
  Name "gt extractseq -keys from large fastafile"
  Keywords "gt_extractseq"
  Test do
    run_test "#{$bin}gt extractseq -o gi-extract.fna.gz -gzip " +
             "-keys #{$gttestdata}gi-queries/gi-queries.txt " +
             "#{$testdata}at1MB"
  end
  Name "gt extractseq -keys from fastaindex"
  Keywords "gt_extractseq"
  Test do
    run_test("#{$bin}gt suffixerator -protein -ssp -tis -des -sds -kys " +
             "-db #{$gttestdata}trembl/trembl-section.fsa.gz")
    run("gunzip -c #{$gttestdata}trembl/trembl-section.fsa.gz")
    run("mv #{last_stdout} trembl-section.fsa")
    run("#{$scriptsdir}/tr2deskeys.rb #{$gttestdata}trembl/trembl-section.fsa.gz")
    run("mv #{last_stdout} trembl-section.keylist")
    run("#{$scriptsdir}/randlines.rb trembl-section.keylist 1000")
    run("mv #{last_stdout} trembl-section.random-keylist")
    run_test("#{$bin}gt extractseq -keys trembl-section.keylist -width 60 " +
             "trembl-section.fsa.gz")
    run("cmp -s #{last_stdout} trembl-section.fsa")
    run_test("#{$bin}gt extractseq -keys #{$testdata}trkeys.txt -width 60 " +
             "trembl-section.fsa.gz")
    run("cmp -s #{last_stdout} #{$testdata}trkeys-result.txt")
    run_test("#{$bin}gt extractseq -keys #{last_stdout} " +
             "trembl-section.fsa.gz",:retval => 1)
    run_test("#{$bin}gt extractseq -keys #{$testdata}trembl-wrongkey.txt " +
             "trembl-section.fsa.gz",:retval => 1)
    run_test("#{$bin}gt extractseq -keys trembl-section.random-keylist -width 60 " +
             "trembl-section.fsa.gz")
    run_test("#{$bin}gt suffixerator -protein -ssp -tis -des -sds -kys sort " +
             "-db #{last_stdout}")
    run("mv #{last_stdout} trembl-section-sorted.fna")
    run_test("#{$bin}gt suffixerator -protein -ssp -tis -des -sds -kys " +
             "-db trembl-section-sorted.fna")
    run("#{$scriptsdir}/tr2deskeys.rb trembl-section-sorted.fna")
    run("mv #{last_stdout} trembl-section-sorted.keylist")
    run_test("#{$bin}gt extractseq -keys trembl-section-sorted.keylist " +
             "-width 60 trembl-section-sorted.fna")
    run("cmp -s #{last_stdout} trembl-section-sorted.fna")
  end
end
