def checkmapped(parts,extraargs,filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run "#{$bin}gt suffixerator -parts #{parts} #{extraargs} -tis -suf -indexname sfx " +
       filearg
  run "#{$bin}gt dev sfxmap sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "cmp -s sfx.prj #{$last_stdout}"
end

1.upto(3) do |parts|
  1.upto(2) do |withsmap|
    if withsmap == 1
      extra="-smap TransProt11"
    else
      extra=""
    end
    Name "gt suffixerator protein #{extra} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checkmapped(parts,"-pl 3 " + extra,["sw100K1.fna","sw100K2.fna"])
    end
  end
end

1.upto(3) do |parts|
  1.upto(2) do |withsmap|
    if withsmap == 1
      extra="-smap TransDNA"
    else
      extra=""
    end
    Name "gt suffixerator dna #{extra} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checkmapped(parts,"-pl 3 " + extra,["Random-Small.fna"])
      checkmapped(parts,"-pl 3 " + extra,["Random.fna"])
      checkmapped(parts,"-pl 3 " + extra,["RandomN.fna"])
      checkmapped(parts,"-pl 3 " + extra,["Random.fna","Atinsert.fna"])
      checkmapped(parts,"-pl 3 " + extra,["RandomN.fna","Random.fna","Atinsert.fna"])
      checkmapped(parts,"-pl 7 " + extra,["RandomN.fna","Random.fna","Atinsert.fna"])
    end
  end
end
