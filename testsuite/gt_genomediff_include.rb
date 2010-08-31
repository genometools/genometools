
fp = File.open("#{$testdata}/genomediff/testsuite", 'r')
smallfiles = fp.readlines
fp.close
fp = File.open("#{$testdata}/genomediff/codelist", 'r')
allfiles = fp.readlines
fp.close
allfiles += smallfiles
kr_testable_files = []
allfiles.each do |code|
  if code.match /[^_]_0+1_/
    kr_testable_files << code.chomp
  end
end

smallfiles.each do |code|
  code.chomp!
  Name "gt genomediff kr2-algorithm traverse"
  Keywords "gt_genomediff kr2-algorithm traverse"
  Test do
    run_test(
      "#{$bin}gt genomediff -traverse -pck #{$testdata}/genomediff/#{code}_idx")
  end
  Name "gt genomediff kr2-algorithm shulen"
  Keywords "gt_genomediff kr2-algorithm shulen"
  Test do
    run_test(
      "#{$bin}gt genomediff -shulen -pck #{$testdata}/genomediff/#{code}_idx")
  end
  Name "gt genomediff kr2-algorithm fullrun"
  Keywords "gt_genomediff kr2-algorithm fullrun"
  Test do
    run_test(
      "#{$bin}gt genomediff -pck #{$testdata}/genomediff/#{code}_idx")
  end
end

kr_testable_files.each do |code|
  Name "gt genomediff kr2-algorith check shulen"
  Keywords "gt_genomediff kr2-algorithm check_shulen"
  Test do
    run_test(
      "#{$bin}gt genomediff -shulen -pck #{$testdata}/genomediff/#{code}_idx",
      :maxtime => 240)
    gt_out = []
    kr_out = []
    numoffiles = 0
    File.open($last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        if numoffiles == 0
          numoffiles = line.match(/^\d+$/)[0].to_i
        else
          gt_out << line.split
        end
      end
      if numoffiles == NIL or
        numoffiles != gt_out.length
        failtest("can't parse output") 
      end
    end
    File.open("#{$testdata}/genomediff/#{code}-kr.out", 'r') do |krfile|
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
    0.upto(numoffiles - 1) do |i|
      1.upto(numoffiles) do |j|
        next if i==j-1
        if gt_out[i][j].to_i - 2 != kr_out[i][j].to_i
          failtest("different results at #{gt_out[i][0]},#{j}")
        end
      end
    end
  end
end

kr_testable_files.each do |code|
  Name "gt genomediff kr2-algorith check kr"
  Keywords "gt_genomediff kr2-algorithm check_kr"
  Test do
    run_test(
      "#{$bin}gt genomediff -pck #{$testdata}/genomediff/#{code}_idx",
      :maxtime => 240)
    gt_out = []
    kr_out = []
    numoffiles = 0
    File.open($last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        if numoffiles == 0
          numoffiles = line.match(/^\d+$/)[0].to_i
        else
          gt_out << line.split
        end
      end
      if numoffiles == NIL or
        numoffiles != gt_out.length
        failtest("can't parse output") 
      end
    end
    File.open("#{$testdata}/genomediff/#{code}-kr.out", 'r') do |krfile|
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
    0.upto(numoffiles - 1) do |i|
      1.upto(numoffiles) do |j|
        next if i==j-1
        next if gt_out[i][j].to_f > 0.3 or kr_out[i][j].to_f > 0.3
        if (gt_out[i][j].to_f - kr_out[i][j].to_f).abs > 0.001
          failtest("different results at #{gt_out[i][0]},#{j}\n"+
                   "#{gt_out[i][j]} #{kr_out[i][j]}")
        end
      end
    end
  end
end
