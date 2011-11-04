def runmerge(queryfile,referencefiles)
  sfxopts="-dna -algbds 3 31 80 -suf -lcp -tis -pl"
  run_test "#{$bin}gt suffixerator #{sfxopts} -indexname all " +
           "-db #{referencefiles.join(" ")}"
  num = 0
  indexlist = Array.new()
  referencefiles.each do |filename|
    run_test "#{$bin}gt suffixerator #{sfxopts} -indexname midx#{num} " +
             "-db #{filename}"
    run_test "#{$bin}gt dev sfxmap -tis -suf -lcp -esa midx#{num}"
    indexlist.push("midx#{num}")
    num+=1
  end
  run_test "#{$bin}gt dev mergeesa -indexname midx-all " +
           "-ii #{indexlist.join(" ")}"
  run "cmp -s midx-all.suf all.suf"
  run "cmp -s midx-all.lcp all.lcp"
  run "cmp -s midx-all.llv all.llv"
  run_test "#{$bin}gt mkfmindex -noindexpos -fmout fm-all " + 
           "-ii #{indexlist.join(" ")}"
  run_test "#{$bin}gt suffixerator -indexname fm-all -plain -des no -ssp no" +
           " -sds no -smap fm-all.al1 -tis -db fm-all.bwt"
  run_test "#{$bin}gt uniquesub -fmi fm-all -query #{queryfile} " +
           "-output sequence querypos -min 10 -max 10"
end

def iterrunmerge(numtoselect)
  reference="#{$testdata}at1MB"
  run "#{$scriptsdir}seqselect.rb #{numtoselect} #{reference}"
  run "#{$scriptsdir}splitmultifasta.rb TMP 0 #{last_stdout}"
  referencefiles = Array.new()
  Dir.new('.').each do |filename|
    if filename.match(/^TMP-/)
      referencefiles.push(filename)
    end
  end
  queryfile="#{$testdata}U89959_genomic.fas"
  runmerge(queryfile,referencefiles)
end

2.upto(10) do |numtoselect|
  Name "gt merge #{numtoselect} enhanced suffix arrays"
  Keywords "gt_mergeesa"
  Test do
    iterrunmerge(numtoselect)
  end
end
