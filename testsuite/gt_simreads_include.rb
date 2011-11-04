def total_length(fas)
  run_test %Q(grep -v ">" #{fas} | tr -d "[:space:]" | wc -c)
  return Integer(IO.read(last_stdout))
end

def nof_seq(fas)
  run_test %Q(grep  "^>" #{fas} | wc -l)
  return Integer(IO.read(last_stdout))
end

def expect(label, value, expected)
  if value != expected
    failtest("#{label}: expected value #{expected}, "+
                       "found #{value}")
  end
end

def expect_range(label, value, min, max)
  if value < min or value > max
    failtest("#{label}: expected value in range #{min}-#{max}, "+
                       "found #{value}")
  end
end

Name "gt simreads -num -len test"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  run_test "#{$bin}gt simreads -num 100 -len 100 "+
                   "-force -o reads #{File.basename(fas)}"

  expect("number of reads", nof_seq("reads"), 100)  
  expect("total length", total_length("reads"), 10000)
end


Name "gt simreads test -minlen/-maxlen"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  run_test "#{$bin}gt simreads -num 100 -minlen 10 -maxlen 100 "+
                   "-force -o reads #{File.basename(fas)}"

  expect("number of reads", nof_seq("reads"), 100)  
  expect_range("total length", total_length("reads"), 1000, 10000)
end

Name "gt simreads -coverage test"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  cov = 5
  len = 100
  run_test "#{$bin}gt simreads -coverage #{cov} -len #{len} "+
                   "-force -o reads #{File.basename(fas)}"
  n = total_length(fas)
  exptlen = n * cov
  mod = exptlen % len
  exptlen += (len-mod) if mod > 0
  expnreads = exptlen / len
  expect("total length", total_length("reads"), exptlen)
  expect("number of reads", nof_seq("reads"), expnreads)
end

Name "gt simreads grep test"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  cov = 3
  len = 100
  n = total_length(fas)

  run_test "#{$bin}gt simreads -coverage #{cov} -len #{len} "+
                   "-force -o reads #{File.basename(fas)}"

  # test that reads originate from target sequence
  # or its reverse complement
  run_test "#{$bin}gt convertseq -fastawidth #{n+1} #{fas}"
  run_test "tail -n 1 #{last_stdout}"
  direct = IO.read(last_stdout)
  run_test "#{$bin}gt convertseq -r -fastawidth #{n+1} #{fas}"
  run_test "tail -n 1 #{last_stdout}"
  rc = IO.read(last_stdout)

  run_test "#{$bin}gt convertseq -fastawidth #{len+1} -force -o reads2 reads"
  
  f = File.open("reads2")
  i = 0
  f.each do |line|
    if !(line =~ /^>/)
      pattern = Regexp.new(Regexp.escape(line.chomp!), Regexp::IGNORECASE)
      if !(direct =~ pattern) && !(rc =~ pattern)
        failtest("read #{i} does not match")
      end
    end
    i+=1
  end  
end

def test_distfile_format(filename)
  f = File.open(filename)
  f.each do |line|
    if (line !~ /^#/ && line.chomp !~ /^\d+ \d+$/)
      failtest("distribution output file format error: #{line}")
    end
  end
  f.close
end

Name "gt simreads -ds test"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  cov = 10
  len = 10
  run_test "#{$bin}gt simreads -coverage #{cov} -len #{len} "+
                  "-ds starts -force -o reads #{File.basename(fas)}"
  test_distfile_format("starts")
end

Name "gt simreads -dl test"
Keywords "gt_simreads"
Test do
  fas = "#{$testdata}U89959_genomic.fas"
  run_test "#{$bin}gt encseq encode #{fas}"

  cov = 10
  minlen = 10
  maxlen = 100
  run_test "#{$bin}gt simreads -coverage #{cov} -minlen #{minlen} "+
                  "-maxlen #{maxlen} -dl lengths -force -o reads #{File.basename(fas)}"
  test_distfile_format("lengths")
end
