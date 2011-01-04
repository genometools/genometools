Name "gt seqencode|seqdecode simple test"
Keywords "gt_seqencode gt_seqdecode"
Test do
  run "#{$bin}gt dev seqencode #{$testdata}foobar.fas"
  run "#{$bin}gt dev seqdecode foobar.fas"
  run "diff #{$last_stdout} #{$testdata}foobar.fas"
end

Name "gt seqencode multiple files without indexname"
Keywords "gt_seqencode"
Test do
  run "#{$bin}gt dev seqencode #{$testdata}foobar.fas"
  run_test "#{$bin}gt dev seqencode #{$testdata}foobar.fas " + \
           "#{$testdata}foobar.fas", :retval => 1
  grep($last_stderr, /if more than one input file is given/)
end

Name "gt seqdecode lossless without ois"
Keywords "gt_seqdecode lossless"
Test do
  run "#{$bin}gt dev seqencode #{$testdata}foobar.fas"
  run_test "#{$bin}gt dev seqdecode -lossless foobar.fas", \
           :retval => 1
  grep($last_stderr, /cannot open file.*ois/)
end

STDREADMODES  = ["fwd", "rev"]
DNAREADMODES  = STDREADMODES + ["cpl", "rcl"]
DNATESTSEQS   = ["#{$testdata}foobar.fas",
                 "#{$testdata}gt_bioseq_succ_3.fas",
                 "#{$testdata}at1MB"]
AATESTSEQS    = ["#{$testdata}trembl-eqlen.faa"]
NUMSAMPLES    = 5

def revcomp(seq)
  comp(seq).reverse
end

def comp(seq)
  seq.tr("aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV","tTgGcCaAnNyYrRkKmMwWsSvVhHdDbB")
end

def mapwildcards(seq)
  seq.downcase.tr("rRyYmMkKwWsSbBdDhHvV", "nNnNnNnNnNnNnNnNnNnN")
end

def getseq(filename, mirrored = false, rm = "fwd")
  sequences = []
  rcseqs = []
  File.open(filename) do |file|
    seqs = file.read.split(">")
    seqs.shift
    seqs.each do |seq|
      seqarr = seq.split(/\n/)
      header = seqarr.shift.chomp
      seq = seqarr.collect{|v| v.chomp}.join('')
      sequences.push(seq)
      rcseqs.push(revcomp(seq))
    end
  end
  ret =  sequences.join("|")
  if mirrored then
    ret += "|" + rcseqs.reverse.join("|")
  end
  case rm
    when "fwd" then
      #pass
    when "rev" then
      ret = ret.reverse
    when "cpl" then
      ret = comp(ret)
    when "rcl" then
      ret = revcomp(ret)
    else
      raise "unknown readmode"
  end
  ret
end

def run_encseq_comparison(filename, mirrored, lossless, readmode, singlechars,
                          numsamples = NUMSAMPLES)
  seq = getseq(filename, mirrored, readmode)
  ranges = []
  numsamples.times do
    len = rand((seq.length)/2)
    start = rand(seq.length-1-len);
    stop = start + len
    ranges.push([start, stop])
  end

  ranges.each do |rng|
    line = "#{$bin}gt dev seqdecode -output concat " + \
           "-range #{rng[0]} #{rng[1]} " + \
           "#{"-lossless" if lossless} " + \
           "-dir #{readmode} #{"-mirrored" if mirrored} " + \
           "#{"-singlechars" if singlechars} #{filename.split('/').last}"
    if mirrored and AATESTSEQS.include?(filename)
      # -mirroring should fail on proteins
      run_test(line, :retval => 1)
    else
      run_test line
      File.open("seqout", "w+") do |f|
        outseq = seq[rng[0]..rng[1]]
        if DNATESTSEQS.include?(filename) and !lossless then
          outseq = mapwildcards(outseq)
        end
        f.write(outseq)
        f.write("\n")
      end
      run "diff seqout #{$last_stdout}"
    end
  end
end

def testformirrored(s, readmode)
  [false, true].each do |lossless|
    [false, true].each do |mirrored|
      [false, true].each do |singlechars|
        Name "gt seqdecode #{s.split('/').last} cc " + \
             "#{"m " if mirrored}#{"s " if singlechars}#{"l " if lossless}#{readmode}"
        Keywords "gt_seqdecode #{" mirroring" if mirrored}#{" lossless" if lossless}"
        Test do
          run "#{$bin}gt dev seqencode -des -ssp -sds " + \
              "#{"-lossless" if lossless} " + \
              "#{s}"
          run_encseq_comparison(s, mirrored, lossless, readmode, singlechars)
        end

        Name "gt seqdecode #{s.split('/').last} cc " + \
             "#{"m " if mirrored}#{"s " if singlechars}#{"l " if lossless}#{readmode} " + \
             "whole seq"
        Keywords "gt_seqdecode#{" mirroring" if mirrored}#{" lossless" if lossless}"
        Test do
          run_test "#{$bin}gt dev seqencode -des -ssp -sds " + \
                   "#{"-lossless" if lossless} " + \
                   "#{s}"
          seq = getseq(s, mirrored, readmode)
          line = "#{$bin}gt dev seqdecode -output concat -dir #{readmode} " + \
                 "#{"-mirrored" if mirrored} " + \
                 "#{"-lossless" if lossless} " + \
                 "#{"-singlechars" if singlechars} ./#{s.split('/').last}"
          if mirrored and AATESTSEQS.include?(s)
            # -mirroring should fail on proteins
            run_test(line, :retval => 1)
          else
            run_test line
            File.open("seqout", "w+") do |f|
              if DNATESTSEQS.include?(s) and !lossless then
                seq = mapwildcards(seq)
              end
              f.write(seq)
              f.write("\n")
            end
            run "diff seqout #{$last_stdout}"
          end
        end
      end
    end
  end
end

DNATESTSEQS.each do |s|
  DNAREADMODES.each do |readmode|
    testformirrored(s, readmode)
  end
end

AATESTSEQS.each do |s|
  STDREADMODES.each do |readmode|
    testformirrored(s, readmode)
  end
end

