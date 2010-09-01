
kr_testable_files = []
smallfilecodes = []
bigfilecodes = []

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

fp = File.open("#{$testdata}/genomediff/testsuite", 'r')
smallfilecodes = fp.readlines
fp.close

smallfilecodes.collect! do |filecode|
  "#{$testdata}/genomediff/#{filecode}"
end

if $gttestdata
  fp = File.open("#{$gttestdata}/genomediff/codelist", 'r')
  bigfilecodes = fp.readlines
  fp.close

  bigfilecodes.collect! do |filecode|
    "#{$gttestdata}/genomediff/#{filecode}"
  end
end

allfilecodes = smallfilecodes + bigfilecodes
allfilecodes.each do |code|
  if code.match /[^_]_0+1_/
    kr_testable_files << code.chomp
  end
end

def test_pck(files, param)
  run_test("#{$bin}gt       " +
           "packedindex mkindex " +
           "-db #{files}        " +
           "-dna                " +
           "-dir rev            " +
           "-ssp                " +
           "-dc 64              " +
           "-bsize 8            " +
           "-sprank             " +
           "-pl                 " +
           "-indexname pck")
  run_test(
    "#{$bin}gt genomediff #{param} -pck pck",
    :maxtime => 240)
end

def test_esa(files)
    run_test "#{$bin}gt suffixerator -db #{files} -indexname esa " + 
             "-dna -suf -tis -lcp"
    run_test "#{$bin}gt shulengthdist -ii esa"
end

def compare_2d_result(matrix1, matrix2)
  0.upto(matrix1.length - 1) do |i_idx|
    1.upto(matrix2.length) do |j_idx|
      if i_idx!=j_idx - 1 and 
         matrix1[i_idx][j_idx] != matrix2[i_idx][j_idx]
        return [i_idx,j_idx]
      end
    end
  end
  return false
end

allfilecodes.each do |code|
  code.chomp!
  Name "gt genomediff kr2-algorithm traverse #{code}"
  Keywords "gt_genomediff kr2-algorithm traverse"
  Test do
    test_pck("#{code}*plus*.fas", "-traverse")
  end
  Name "gt genomediff kr2-algorithm shulen #{code}"
  Keywords "gt_genomediff kr2-algorithm shulen"
  Test do
    test_pck("#{code}*plus*.fas", "-shulen")
  end
  Name "gt genomediff kr2-algorithm fullrun #{code}"
  Keywords "gt_genomediff kr2-algorithm fullrun"
  Test do
    test_pck("#{code}*plus*.fas", "")
  end
end

kr_testable_files.each do |code|
  Name "gt genomediff kr2-algorith check shulen #{code}"
  Keywords "gt_genomediff kr2-algorithm check_shulen"
  Test do
    pck_out = []
    esa_out =[]
    kr_out = []
    numoffiles = 0

    test_pck("#{code}*plus*.fas", "-shulen")

    File.open($last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        if numoffiles == 0
          numoffiles = line.match(/^\d+$/)[0].to_i
        else
          pck_out << line.split
        end
      end
      if numoffiles == nil or
        numoffiles != pck_out.length
        failtest("can't parse output") 
      end
    end

    test_esa("#{code}*plus*.fas")

    File.open($last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        unless line.match(/^\d+$/)
          esa_out << line.split
        end
      end
      if numoffiles == nil or
        numoffiles != pck_out.length
        failtest("can't parse output or wrong line numbers") 
      end
    end

    File.open("#{code}-kr.out", 'r') do |krfile|
      line = krfile.gets
      line.chomp!
      if numoffiles != line.to_i
        failtest("different num of files")
      end
      while line = krfile.gets do
        line.chomp!
        if line.match /^\d+$/
          break
        end
        kr_out << line.split
      end
    end
    if result = compare_2d_result(pck_out,kr_out)
      failtest("different results pck-kr #{result[0]},#{result[1]}")
    end
    if result = compare_2d_result(pck_out,esa_out)
      failtest("different results pck-esa #{result[0]},#{result[1]}")
    end
    if result = compare_2d_result(esa_out,kr_out)
      failtest("different results esa-kr #{result[0]},#{result[1]}")
    end
  end
end

kr_testable_files.each do |code|
  Name "gt genomediff kr2-algorith check kr #{code}"
  Keywords "gt_genomediff kr2-algorithm check_kr"
  Test do
    pck_out = []
    kr_out = []
    numoffiles = 0

    test_pck("#{code}*plus*.fas", "")

    File.open($last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        if numoffiles == 0
          numoffiles = line.match(/^\d+$/)[0].to_i
        else
          pck_out << line.split
        end
      end
      if numoffiles == nil or
        numoffiles != pck_out.length
        failtest("can't parse output") 
      end
    end
    File.open("#{code}-kr.out", 'r') do |krfile|
      line = krfile.gets
      line.chomp!
      if numoffiles != line.to_i
        failtest("different num of files")
      end
      (numoffiles+1).times do
        krfile.gets
      end
      while line = krfile.gets do
        line.chomp!
        if line.match /^\d+$/
          break
        end
        kr_out << line.split
      end
    end
    0.upto(numoffiles - 1) do |i_idx|
      1.upto(numoffiles) do |j_idx|
        unless i_idx==j_idx-1 or
           pck_out[i_idx][j_idx].to_f > 0.3 or
           kr_out[i_idx][j_idx].to_f > 0.3
          if (pck_out[i_idx][j_idx].to_f -
              kr_out[i_idx][j_idx].to_f).abs > 0.001
            failtest("different results at #{pck_out[i_idx][0]},#{j_idx}\n"+
                     "#{pck_out[i_idx][j_idx]} #{kr_out[i_idx][j_idx]}")
          end
        end
      end
    end
  end
end

def check_shulen_for_list_pairwise(list)
  list.each do |file1|
    list.each do |file2|
      if file1 != file2
        numoffiles = 0
        pck_out = []
        esa_out =[]

        test_pck("#{$testdata}/#{file1} #{$testdata}/#{file2}", "-shulen")

        File.open($last_stdout, 'r') do |outfile|
          while line = outfile.gets do
            line.chomp!
            next if line.match /^#/
            if numoffiles == 0
              numoffiles = line.match(/^\d+$/)[0].to_i
            else
              pck_out << line.split
            end
          end
          if numoffiles == nil or
            numoffiles != pck_out.length
            failtest("can't parse output pck #{file1} #{file2}") 
          end
        end
        test_esa("#{$testdata}/#{file1} #{$testdata}/#{file2}")

        File.open($last_stdout, 'r') do |outfile|
          while line = outfile.gets do
            line.chomp!
            next if line.match /^#/
            unless line.match(/^\d+$/)
              esa_out << line.split
            end
          end
          if numoffiles == NIL or
            numoffiles != pck_out.length
            failtest("can't parse output esa #{file1} #{file2}") 
          end
        end
        if result = compare_2d_result(pck_out,esa_out)
          failtest("different results pck-esa #{result[0]},#{result[1]}")
        end
      end
    end
  end
end

Name "gt genomediff esa pck pairwise smallfiles"
Keywords "gt_genomediff small pairwise check_shulen"
Test do
  check_shulen_for_list_pairwise(allfiles)
end

Name "gt genomediff esa pck pairwise bigfiles"
Keywords "gt_genomediff big pairwise check_shulen"
Test do
  check_shulen_for_list_pairwise(bigfiles)
end

Name "gt genomediff small multi pck esa"
Keywords "gt_genomediff small check_shulen"
Test do
  realfiles = ""
  allfiles.each do |file|
    realfiles += "#{$testdata}/"+ file + " "
  end
  numoffiles = 0
  pck_out = []
  esa_out =[]

  test_pck(realfiles, "-shulen")

  File.open($last_stdout, 'r') do |outfile|
    while line = outfile.gets do
      line.chomp!
      next if line.match /^#/
      if numoffiles == 0
        numoffiles = line.match(/^\d+$/)[0].to_i
      else
        pck_out << line.split
      end
    end
    if numoffiles == nil or
      numoffiles != pck_out.length
      failtest("can't parse output") 
    end
  end
  test_esa(realfiles)

  File.open($last_stdout, 'r') do |outfile|
    while line = outfile.gets do
      line.chomp!
      next if line.match /^#/
      unless line.match(/^\d+$/)
        esa_out << line.split
      end
    end
    if numoffiles == NIL or
      numoffiles != pck_out.length
      failtest("can't parse output or wrong line numbers") 
    end
  end
  if result = compare_2d_result(pck_out,esa_out)
    failtest("different results pck-esa #{result[0]},#{result[1]}")
  end
end
