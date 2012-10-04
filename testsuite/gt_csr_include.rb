require 'fileutils'
require 'tempfile'

hcr_testfiles = [
  "csr_testcase.fastq",
  "description_test.fastq",
  "description_test2.fastq"
]

Name "gt hcr reads"
Keywords "gt_csr hcr_nodesc"
Test do
  hcr_testfiles.each do |file|
    run_test "#$bin/gt compreads compress -files #$testdata/#{file} -name test"
    run_test "#$bin/gt compreads decompress -file test"
    `grep -v @ #$testdata/#{file} > original`
    `grep -v @ test.fastq > test_out`
    run_test "diff test_out original"
  end
end

Name "gt hcr reads allfiles"
Keywords "gt_csr hcr_nodesc"
Test do
  files = hcr_testfiles.collect{|file| "#$testdata/" + file}
  run_test "#$bin/gt compreads compress -files #{files.join(' ')} -name test"
  run_test "#$bin/gt compreads decompress -file test"
  `grep -h -v @ #{files.join(' ')} > original`
  `grep -v @ test.fastq > test_out`
  run_test "diff test_out original"
end

Name "gt hcr reads and description"
Keywords "gt_csr hcr_desc"
Test do
  hcr_testfiles.each do |file|
    run_test "#$bin/gt compreads compress -descs" +
             " -files #$testdata/#{file} -name test"
    run_test "#$bin/gt compreads decompress -file test"
    `grep -v @ #$testdata/#{file} > original`
    `grep -v @ test.fastq > test_out`
    run_test "diff test_out original"
    run_test "#$bin/gt compreads decompress -descs -file test"
    run_test "diff test.fastq #$testdata/#{file}"
  end
end

Name "gt hcr reads and description allfiles"
Keywords "gt_csr hcr_desc"
Test do
  files = hcr_testfiles.collect{|file| "#$testdata/" + file}
  run_test "#$bin/gt compreads compress -descs" +
           " -files #{files.join(' ')} -name test"
  run_test "#$bin/gt compreads decompress -descs -file test"
  `cat #{files.join(' ')} > original`
  run_test "diff test.fastq original"
end

Name "gt hcr decompress benchmark"
Keywords "gt_csr hcr benchmark"
Test do
  hcr_testfiles.each do |file|
    run_test "#$bin/gt compreads compress -descs" +
             " -files #$testdata/#{file} -name test_#{file}"
    run_test "#$bin/gt compreads decompress -descs -benchmark 10000" +
             " -file test_#{file}", :maxtime => 300
  end
end

rcr_testfiles = {
  "rcr_testreads_on_seq.bam" => "rcr_testseq.fa",
  "example_1.sorted.bam" => "example_1.fa"
}

Name "gt rcr reads noqual"
Keywords "gt_csr rcr"
Test do
  rcr_testfiles.keys do |file|
    run_test "#$bin/gt encseq encode -dna" +
             " -indexname ./#{rcr_testfiles[file]}" +
             " #$testdata/#{rcr_testfiles[file]}"
    run_test "#$bin/gt compreads refcompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -bam #$testdata/#{file}" +
             " -name #{file}"
    run_test "#$bin/gt compreads refdecompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -rcr ./#{file}"
  end
end

Name "gt rcr reads qual"
Keywords "gt_csr rcr"
Test do
  rcr_testfiles.keys do |file|
    run_test "#$bin/gt encseq encode -dna" +
             " -indexname ./#{rcr_testfiles[file]}" +
             " #$testdata/#{rcr_testfiles[file]}"
    run_test "#$bin/gt compreads refcompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -bam #$testdata/#{file}" +
             " -mquals -quals" +
             " -name #{file}"
    run_test "#$bin/gt compreads refdecompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -rcr ./#{file}"
  end
end

Name "gt rcr reads variant qual"
Keywords "gt_csr rcr"
Test do
  rcr_testfiles.keys do |file|
    run_test "#$bin/gt encseq encode -dna" +
             " -indexname ./#{rcr_testfiles[file]}" +
             " #$testdata/#{rcr_testfiles[file]}"
    run_test "#$bin/gt compreads refcompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -bam #$testdata/#{file}" +
             " -mquals -vquals" +
             " -name #{file}"
    run_test "#$bin/gt compreads refdecompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -rcr ./#{file}"
  end
end

Name "gt rcr reads variant qual, descriptions"
Keywords "gt_csr rcr"
Test do
  rcr_testfiles.keys do |file|
    run_test "#$bin/gt encseq encode -dna" +
             " -indexname ./#{rcr_testfiles[file]}" +
             " #$testdata/#{rcr_testfiles[file]}"
    run_test "#$bin/gt compreads refcompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -bam #$testdata/#{file}" +
             " -mquals -vquals -descs" +
             " -name #{file}"
    run_test "#$bin/gt compreads refdecompress" +
             " -ref ./#{rcr_testfiles[file]}" +
             " -rcr ./#{file}"
             " -qnames"
  end
end
