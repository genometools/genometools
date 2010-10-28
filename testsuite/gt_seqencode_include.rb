Name "gt seqencode|seqdecode test"
Keywords "gt_seqencode gt_seqdecode"
Test do
  run "#{$bin}gt dev seqencode #{$testdata}foobar.fas"
  run "#{$bin}gt dev seqdecode #{$testdata}foobar.fas"
  run "diff #{$last_stdout} #{$testdata}foobar.fas"
end

STDREADMODES  = ["fwd", "rev"]
DNAREADMODES  = STDREADMODES + ["cpl", "rcl"]
DNATESTSEQS   = ["#{$testdata}foobar.fas",
                 "#{$testdata}gt_bioseq_succ_3.fas",
                 "#{$testdata}at1MB"]
AATESTSEQS    = ["#{$testdata}trembl-eqlen.faa"]
NUMSAMPLES    = 10

def revcomp(seq)
  comp(seq).reverse
end

def comp(seq)
  seq.tr("aAcCgGtTnN","tTgGcCaAnN")
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

def run_encseq_comparison(filename, mirrored, readmode, singlechars,
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
           "-dir #{readmode} #{"-mirrored" if mirrored} " + \
           "#{"-singlechars" if singlechars} #{filename}"
    if mirrored and AATESTSEQS.include?(filename)
      # -mirroring should fail on proteins
      run_test(line, :retval => 1)
    else
      run_test line
      File.open("seqout", "w+") do |f|
        outseq = seq[rng[0]..rng[1]]
        if DNATESTSEQS.include?(filename)
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
  [false, true].each do |mirrored|
    [false, true].each do |singlechars|
      Name "gt seqdecode #{s.split('/').last} cc " + \
           "#{"mirr " if mirrored}#{"sngl " if singlechars}#{readmode}"
      Keywords "gt_seqdecode #{" mirroring" if mirrored}"
      Test do
        run "#{$bin}gt dev seqencode #{s}"
        run_encseq_comparison(s, mirrored, readmode, singlechars)
      end

      Name "gt seqdecode #{s.split('/').last} cc " + \
           "#{"mirr " if mirrored}#{"sngl " if singlechars}#{readmode} " + \
           "whole seq"
      Keywords "gt_seqdecode#{" mirroring" if mirrored}"
      Test do
        run_test "#{$bin}gt dev seqencode #{s}"
        seq = getseq(s, mirrored, readmode)
        line = "#{$bin}gt dev seqdecode -output concat -dir #{readmode} " + \
               "#{"-mirrored" if mirrored} " + \
               "#{"-singlechars" if singlechars} #{s}"
        if mirrored and AATESTSEQS.include?(s)
          # -mirroring should fail on proteins
          run_test(line, :retval => 1)
        else
          run_test line
          File.open("seqout", "w+") do |f|
            if DNATESTSEQS.include?(s)
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

