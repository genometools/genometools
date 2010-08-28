allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "TTT-small.fna",
            "trna_glutamine.fna",
            "Atinsert.fna",
            "U89959_genomic.fas"]

def checkshulengthdist(file1,file2)
  run_test "#{$bin}gt suffixerator -db #{file1} -indexname index1 " + 
           "-dna -suf -tis"
  run_test "#{$bin}gt shulengthdist -ii index1 -q #{file2}"
  run "mv #{$last_stdout} result.pairwise"
  run_test "#{$bin}gt suffixerator -db #{file2} -indexname index2 " + 
           "-dna -suf -tis"
  run_test "#{$bin}gt shulengthdist -ii index2 -q #{file1}"
  run "cat #{$last_stdout} >> result.pairwise"
  run_test "#{$bin}gt suffixerator -db #{file1} #{file2} -indexname both " + 
           "-dna -suf -tis -lcp"
  run_test "#{$bin}gt shulengthdist -ii both"
  run "mv #{$last_stdout} result.multi"
  run "diff result.pairwise result.multi"
end

allfiles.each do |file1|
  allfiles.each do |file2|
    if file1 != file2
      Name "gt shulengthdist #{file1} #{file2}"
      Keywords "gt_shulengthdist small"
      Test do
        checkshulengthdist("#{$testdata}/#{file1}",
                           "#{$testdata}/#{file2}")
      end
    end
  end
end
