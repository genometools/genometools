require 'fileutils'
require 'tempfile'
require '../scripts/genomediff.rb'

Name "gt genomediff pck <esa> failure"
Keywords "gt_genomediff pck failure"
Test do
  run_test "#{$bin}gt suffixerator -suf -lcp -indexname esa " +
           "-db #{$testdata}Atinsert.fna"
  run_test "#{$bin}gt genomediff -indextype pck esa",:retval=>1
end

Name "gt genomediff esa <pck> failure"
Keywords "gt_genomediff esa failure"
Test do
  run_test "#{$bin}gt packedindex mkindex -indexname pck " +
           "-db #{$testdata}Atinsert.fna"
  run_test "#{$bin}gt genomediff -indextype esa pck",:retval=>1
end

def reverse_and_concat(file)
  tf = Tempfile.new("gt_rev")
  tmpfile = tf.path
  tf.close
  dirname = File.dirname(file)
  `#{$bin}gt                \
   convertseq -o #{tmpfile} \
   -force -r #{file}`
  basename = File.basename(file, ".*")
  newfile = basename + "_plus_rev.fas"
  FileUtils.cp(file, newfile)
  file = newfile
  File.open(file, 'a') {|fp|
    File.open(tmpfile, 'r') {|tf|
      while line = tf.gets
        fp.puts line
      end
    }
  }
  failtest unless File.exist?(newfile)
  return newfile
end

kr_testable_files = []
smallfilecodes = []
bigfilecodes = []

allfiles = ["Duplicate.fna",
            "Random159.fna",
            "Random160.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

bigfiles = ["U89959_genomic.fas",
            "Random.fna",
            "Atinsert.fna"]

# do not run the big tests with valgrind
if $gttestdata and not $arguments["memcheck"]
  bigfiles.push "at1MB"
end

fp = File.open("#{$testdata}genomediff/testsuite", 'r')
smallfilecodes = fp.readlines
fp.close

smallfilecodes.collect! do |filecode|
  "#{$testdata}genomediff/#{filecode}".chomp
end

# do not run the big tests with valgrind
if $gttestdata and not $arguments["memcheck"]
  fp = File.open("#{$gttestdata}genomediff/codelist", 'r')
  bigfilecodes = fp.readlines
  fp.close

  bigfilecodes.collect! do |filecode|
    "#{$gttestdata}genomediff/#{filecode}".chomp
  end
end

allfilecodes = smallfilecodes + bigfilecodes
allfilecodes.each do |code|
  if code.match /[^_]_0+1_/
    kr_testable_files << code.chomp
  end
end

def test_pck(files, param, idxparam)
  run_test("#{$bin}gt       " +
           "packedindex mkindex " +
           "-db #{files} " +
           "-dna " +
           "-dir rev " +
           "-ssp " +
           "-dc 64 " +
           "-bsize 8 " +
           "-sprank " +
           "-pl " +
           "#{idxparam} " +
           "-indexname pck")
  run_test(
    "#{$bin}gt genomediff #{param} -indextype pck pck",
    :maxtime => 720)
end

def test_esa(files, param, idxparam)
  run_test "#{$bin}gt suffixerator -db #{files} -indexname esa " +
           "-dna -suf -tis -lcp -ssp #{idxparam}"
  run_test "#{$bin}gt genomediff #{param} -indextype esa esa"
end

def test_esq(files, param)
  run_test "#$bin/gt genomediff #{param} #{files}"
end

def compare_2d_result(matrix1, matrix2)
  0.upto(matrix1.length - 1) do |i_idx|
    1.upto(matrix2.length) do |j_idx|
      if i_idx!=j_idx - 1 and
         matrix1[i_idx][j_idx] != matrix2[i_idx][j_idx]
        return matrix1[i_idx][j_idx], matrix2[i_idx][j_idx]
      end
    end
  end
  return false
end

Name "gt genomediff pck testset"
Keywords "gt_genomediff pck"
Test do
  allfilecodes.each do |code|
      test_pck("#{code}*.fas", "", "")
      test_pck("#{code}*.fas", "", "-mirrored")
  end
end
Name "gt genomediff esa testset"
Keywords "gt_genomediff esa"
Test do
  allfilecodes.each do |code|
    test_esa("#{code}*.fas", "", "")
    test_esa("#{code}*.fas", "", "-mirrored")
  end
end

Name "gt genomediff esq testset"
Keywords "gt_genomediff esq"
Test do
  allfilecodes.each do |code|
    test_esq("#{code}*.fas", "-indexname esq")
    test_esq("#{code}*.fas", "-indexname esq -mirrored")
    test_esq("esq", "")
    test_esq("esq", "-mirrored")
  end
end

code = smallfilecodes[0]

Name "gt genomediff compare kr"
Keywords "gt_genomediff esa pck check_kr simulated_data"
Test do
  kr_testable_files.each do |code|
    pck_out = []
    esa_out = []
    esq_out = []
    kr_out = []
    numoffiles = 0
    files = " "
    Dir.glob("#{code}*.fas").sort.each do |file|
      files += reverse_and_concat(file) + ' '
    end

    test_pck("#{files}", "", "")

    File.open(last_stdout, 'r') do |outfile|
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

    test_esa("#{files}", "", "")

    File.open(last_stdout, 'r') do |outfile|
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

    test_esq("#{files}", "-indexname esq")

    File.open(last_stdout, 'r') do |outfile|
      while line = outfile.gets do
        line.chomp!
        next if line.match /^#/
        unless line.match(/^\d+$/)
          esq_out << line.split
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
    if result = compare_2d_result(pck_out,esa_out)
      failtest("different results pck-esa #{result[0]},#{result[1]}")
    end
    if result = compare_2d_result(pck_out,esq_out)
      failtest("different results pck-esa #{result[0]},#{result[1]}")
    end
    0.upto(numoffiles - 1) do |i_idx|
      1.upto(numoffiles) do |j_idx|
        unless i_idx==j_idx-1 or
           pck_out[i_idx][j_idx].to_f > 0.3 or
           kr_out[i_idx][j_idx].to_f > 0.3
          if (pck_out[i_idx][j_idx].to_f -
              kr_out[i_idx][j_idx].to_f).abs > 0.001
            failtest("different results pck-kr #{pck_out[i_idx][j_idx]},"+
                     "#{kr_out[i_idx][j_idx]}\n"+
                     "#{pck_out[i_idx][j_idx]} #{kr_out[i_idx][j_idx]}")
          end
        end
      end
    end
  end
end

def check_shulen_for_list_pairwise(list)
  Name "gt genomediff pairwise test"
  Keywords "gt_genomediff pairwise esa pck check_shulen"
  Test do
    start_idx = 1
    list.each do |file1|
      list[start_idx..-1].each do |file2|
        if file1 != file2
          numoffiles = 0
          pck_out = []
          esa_out = []

          test_pck("#{$testdata}#{file1} #{$testdata}#{file2}", "", "")

          File.open(last_stdout, 'r') do |outfile|
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
          test_esa("#{$testdata}#{file1} #{$testdata}#{file2}",
                    "", "")

          File.open(last_stdout, 'r') do |outfile|
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
      start_idx += 1
    end
  end
end

check_shulen_for_list_pairwise(allfiles)

check_shulen_for_list_pairwise(bigfiles)

Name "gt genomediff smallfiles all at once"
Keywords "gt_genomediff esa pck small check_shulen"
Test do
  realfiles = ""
  allfiles.each do |file|
    realfiles += "#{$testdata}"+ file + " "
  end
  numoffiles = 0
  pck_out = []
  esa_out =[]

  test_pck(realfiles, "", "")

  File.open(last_stdout, 'r') do |outfile|
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
  test_esa(realfiles, "", "")

  File.open(last_stdout, 'r') do |outfile|
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
