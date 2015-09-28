require "fileutils"

files_dna = ["#{$testdata}Arabidopsis-C99826.fna",
             "#{$testdata}Atinsert.fna",
             "#{$testdata}Duplicate.fna",
             "#{$testdata}example_1.fa",
             "#{$testdata}Random.fna",
             "#{$testdata}RandomN.fna",
             "#{$testdata}Reads1.fna",
             "#{$testdata}Reads2.fna",
             "#{$testdata}Reads3.fna",
             "#{$testdata}Scaffold_102.fa"]

files_protein = ["#{$testdata}nGASP/protein_10.fas",
                 "#{$testdata}nGASP/protein_100_with_stop.fas",
                 "#{$testdata}nGASP/protein_long.fas"]

files_kmersize = ["#{$testdata}TTT-small.fna"]

file_big = []
if $gttestdata
  files_big = ["#{$testdata}condenseq/varlen_0.01_50.fas",
               "#{$testdata}condenseq/varlen_0.01_200.fas"]
end

1.step(7, 3) do |i|
  10.step(110, 50) do |j|
    Name "gt kmer_database dna #{i} #{j}"
    Keywords "gt_kmer_database dna"
    Test do
      #DNA test
      files_dna.each do |file_name|
        FileUtils.copy(file_name, ".")
        run_test "#{$bin}gt encseq encode #{File.basename(file_name)}"
        run_test "#{$bin}gt dev kmer_database -kmersize #{i} -bsize #{j} " \
          "#{File.basename(file_name)}", :maxtime => 300
      end
    end
  end
end
1.step(3, 2) do |i|
  [50, 100, 500].each do |j|
    Name "gt kmer_database prot #{i} #{j}"
    Keywords "gt_kmer_database prot"
    Test do
      #Protein test, kmersize needs to be smaller
      files_protein.each do |file_name|
        FileUtils.copy(file_name, ".")
        run_test "#{$bin}gt encseq encode #{File.basename(file_name)}"
        run_test "#{$bin}gt dev kmer_database -kmersize #{i} -bsize #{j} " \
          "#{File.basename(file_name)}", :maxtime => 300
      end
    end
  end
end

Name "gt kmer_database k too big"
Keywords "gt_kmer_database kmersize fail"
Test do
  FileUtils.copy(files_kmersize[0], ".")
  run_test "#{$bin}gt encseq encode #{File.basename(files_kmersize[0])}"
  run_test "#{$bin}gt dev kmer_database -kmersize 10 " \
    "#{File.basename(files_kmersize[0])}", :retval => 1
  grep last_stderr, "Input is too short for used kmersize. File length:" \
    " 7 kmersize: 10"
end

if $gttestdata
  Name "gt kmer_database large files"
  Keywords "gt_kmer_database large"
  Test do
    files_big.each do |file_name|
      FileUtils.copy(file_name, ".")
      run_test "#{$bin}gt encseq encode #{File.basename(file_name)}"
      run_test "#{$bin}gt dev kmer_database -merge_only -bsize 1000 " \
        " #{File.basename(file_name)}"
    end
  end
end
