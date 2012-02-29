require 'multidimarray.rb'

allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

bigfiles = ["at1MB",
            "U89959_genomic.fas",
            "Atinsert.fna"]

def checkshulengthdistforlist(filelist)
  numofdbfiles = filelist.length
  realfilelist = []
  filelist.each do |filename|
    realfilelist.push("#{$testdata}#{filename}")
  end
  run_test "#{$bin}gt suffixerator -db #{realfilelist.join(' ')} " +
           "-indexname all -dna -suf -lcp"
  run_test "#{$bin}gt genomediff -shulen -esa all"
  run "mv #{last_stdout} shulen.gd"
  run_test "#{$bin}gt suffixerator -db #{realfilelist.join(' ')} " +
           "-indexname all -dna -genomediff -lcp"
  run "mv #{last_stdout} shulen.sfx"
  run "diff shulen.gd  shulen.sfx"
end

allfiles.each do |file1|
  allfiles.each do |file2|
    if file1 != file2
      Name "gt shulengthdist #{file1} #{file2}"
      Keywords "gt_shulengthdist small"
      Test do
        checkshulengthdistforlist([file1,file2])
      end
    end
  end
end

Name "gt shulengthdist allfiles"
Keywords "gt_shulengthdist all"
Test do
  checkshulengthdistforlist(allfiles)
end

Name "gt shulengthdist bigfiles"
Keywords "gt_shulengthdist big"
Test do
  checkshulengthdistforlist(bigfiles)
end
