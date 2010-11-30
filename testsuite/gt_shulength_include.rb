=begin
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

def readnumfromfile(filename)
  f = File.new(filename)
  f.each_line do |line|
    ns = line.match(/(^[0-9]+$)/)
    if ns
      yield ns[1].to_i
    end
  end
  f.close_read
end

def read3numfromfile(filename)
  f = File.new(filename)
  f.each_line do |line|
    ns = line.match(/^([0-9]+)\s+([0-9]+)\s+([0-9]+)$/)
    if ns
      yield ns[1].to_i,ns[2].to_i,ns[3].to_i
    end
  end
end

def compueshulengthdistforpair(file1,file2,indexname)
  run_test "#{$bin}gt suffixerator -db #{file1} -indexname #{indexname} " + 
           "-dna -suf -tis"
  run_test "#{$bin}gt shulengthdist -ii #{indexname} -q #{file2}"
  readnumfromfile($last_stdout) do |n|
    return n
  end
end

def checkshulengthdistforpair(file1,file2)
  dist1 = compueshulengthdistforpair(file1,file2,"index1")
  dist2 = compueshulengthdistforpair(file2,file1,"index2")
  run_test "#{$bin}gt suffixerator -db #{file1} #{file2} -indexname both " + 
           "-dna -suf -tis -lcp"
  run_test "#{$bin}gt shulengthdist -ii both"
  list = []
  read3numfromfile($last_stdout) do |idx1,idx2,n|
    list.push(n)
  end
  if [dist1,dist2] != list
    STDERR.puts "dist1=#{dist1},dist2=#{dist2} != #{list}=list"
    exit 1
  end
end

def checkshulengthdistforlist(filelist)
  numofdbfiles = filelist.length
  realfilelist = []
  filelist.each do |filename|
    realfilelist.push("#{$testdata}#{filename}")
  end
  pairwiseshulengthmatrix = MultiDimensionalArray.new(numofdbfiles,numofdbfiles)
  0.upto(numofdbfiles-1) do |idx|
    pairwiseshulengthmatrix[idx,idx] = 0
  end
  realfilelist.each_with_index do |file1,idx1|
    realfilelist.each_with_index do |file2,idx2|
      if file1 != file2
        dist1 = compueshulengthdistforpair(file1,file2,"index.#{idx1}.#{idx2}")
        pairwiseshulengthmatrix[idx1,idx2] = dist1
        dist2 = compueshulengthdistforpair(file2,file1,"index.#{idx2}.#{idx1}")
        pairwiseshulengthmatrix[idx2,idx1] = dist2
      end
    end
  end
  multishulengthmatrix = MultiDimensionalArray.new(numofdbfiles,numofdbfiles)
  0.upto(numofdbfiles-1) do |idx|
    multishulengthmatrix[idx,idx] = 0
  end
  run_test "#{$bin}gt suffixerator -db #{realfilelist.join(' ')} -indexname all " + 
           "-dna -suf -tis -lcp"
  run_test "#{$bin}gt shulengthdist -ii all"
  read3numfromfile($last_stdout) do |idx1,idx2,n|
    multishulengthmatrix[idx1,idx2] = n
  end
  0.upto(numofdbfiles-1) do |idx1|
    0.upto(numofdbfiles-1) do |idx2|
      pairvalue = pairwiseshulengthmatrix[idx1,idx2]
      multivalue = multishulengthmatrix[idx1,idx2]
      if pairvalue != multivalue
	STDERR.puts "#{idx1},#{idx2}: pairwiseshulength=#{pairvalue} !=" +
                    " #{multivalue} = multishulengthvalue"
        exit 1
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

allfiles.each do |file1|
  allfiles.each do |file2|
    if file1 != file2
      Name "gt shulengthdist #{file1} #{file2}"
      Keywords "gt_shulengthdist small"
      Test do
        checkshulengthdistforpair("#{$testdata}#{file1}",
                                  "#{$testdata}#{file2}")
      end
    end
  end
end
=end
