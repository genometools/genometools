require "set"

Match = Struct.new("Match",:len1,:seq1,:start1,:len2,:seq2,:start2,:score,
                            :distance,:identity,:other,:result1,:result2)

def match_to_s(m)
  return [m.len1,m.seq1,m.start1,m.len2,m.seq2,m.start2,m.score,m.distance].
                      join(" ") + sprintf(" %.2f",m.identity)
end

def prefixlength_get(prjfile)
  prefixlength = nil
  File.open(prjfile).each_line do |line|
    if m = line.match(/^prefixlength=(\d+)$/)
      prefixlength = m[1].to_i
    end
  end
  if prefixlength.nil?
    STDERR.puts "#{prjfile} does not contain definition of prefixlength"
    exit 1
  end
  return prefixlength
end

def createrepfindcall(indexname,seedlength,extend_opt,optionlist,
                      emptyenv,
                      verbose = true)
  verboseoption = if verbose then "-v" else "" end
  repfindcall = "bin/gt repfind -scan #{verboseoption} -minid 80 " +
                "-seedlength #{seedlength} -#{extend_opt} -ii #{indexname} " + 
                optionlist.join(" ")
  if emptyenv 
    return "env -i #{repfindcall}" 
  else 
    return repfindcall 
  end
end

def match_new_without_result(a)
  return Match.new(a[0].to_i,a[1].to_i,a[2].to_i,a[4].to_i,a[5].to_i,a[6].to_i,
                   a[7].to_i,a[8].to_i,a[9].to_f,nil,nil,nil)
end

def makeseedhash(indexname,seedlength,extend_opt,optionlist,gencall = false)
  seedhash = Hash.new()
  key = nil
  if seedlength == 0
    seedlength = 2 * prefixlength_get("#{indexname}.prj")
  end
  repfindcall = createrepfindcall(indexname,seedlength,extend_opt,optionlist,
                                  true,if gencall then false else true end)
  if gencall
    return repfindcall
  end
  STDERR.puts "\# #{repfindcall}"
  IO.popen(repfindcall.split(/\s/)).each_line do |line|
    if line.match(/# seed:/)
      key = line.chomp.gsub(/# seed:\s+/,"")
    elsif line.match(/^#/)
      next
    elsif line.match(/^\d/)
      a = line.split(/\s/)
      if key.nil?
        STDERR.puts "#{$0}: expect that key is defined"
        exit 1
      end
      if seedhash.has_key?(key)
        STDERR.puts "#{$0}: key #{key} already occurs"
        exit 1
      end
      seedhash[key] = match_new_without_result(a)
      key = nil
    end
  end
  if "#{$?}" != "" and not "#{$?}".match(/exit 0$/)
    STDERR.puts "FAILURE: #{repfindcall}: \"#{$?}\""
    exit 1
  end
  return seedhash
end

def inputseedhash(matchfile)
  seedhash = Hash.new()
  seed = nil
  File.open(matchfile,"r").each_line do |line|
    if line.match(/# seed:/)
      seed = line.chomp.gsub(/# seed:\s+/,"")
    elsif line.match(/^#/)
      next
    elsif line.match(/^\d/)
      a = line.split(/\s/)
      if seed.nil?
        STDERR.puts "#{$0}: expect that seed is defined"
        exit 1
      end
      if seedhash.has_key?(seed)
        STDERR.puts "#{$0}: seed #{seed} already occurs"
        exit 1
      end
      seedhash[seed] = match_new_without_result(a)
      seed = nil
    end
  end
  return seedhash
end

def is_contained(start1,end1,start2,end2)
  if start2 <= start1 and end1 <= end2
    return true
  else
    return false
  end
end

def overlap(start1,end1,start2,end2)
  if start1 <= start2
    if start2 <= end1
      return end1 - start2 + 1
    else
      return 0
    end
  elsif start1 <= end2
    return end2 - start1 + 1
  else
    return 0
  end
end

Results = Struct.new("Result",:contained,:containing,:ovperc)

def evaluate_itv(start1,len1,start2,len2)
  end1 = start1 + len1 - 1
  end2 = start2 + len2 - 1
  contained = false
  containing = false
  if is_contained(start1,end1,start2,end2)
    contained = true
  end
  if is_contained(start2,end2,start1,end1)
    containing = true
  end
  ovperc = nil
  if not contained and not containing
    ov = overlap(start1,end1,start2,end2)
    ovperc = 100.0 * ov.to_f/[len1,len2].max.to_f
  end
  return Results.new(contained,containing,ovperc)
end

def addpercentage(h,val)
  sval = sprintf("%.1f",val)
  if h.has_key?(sval)
    h[sval] += 1
  else
    h[sval] = 1
  end
end

def seedhash2seqnum_pairs(seedhash)
  return Set.new(seedhash.values.map {|m| [m.seq1,m.seq2]})
end

def showcomment(size,sum_size,comment)
  if sum_size > 0
    perc = 100.0 * size.to_f/sum_size.to_f
    STDERR.printf("# %d (%.0f%%) of %d sequence pairs %s\n",size,perc,
                                                  sum_size,comment)
    return perc.to_i
  else
    return 0
  end
end

def calcdifference(seqnumpair_set1,seqnumpair_set2)
  sum_size = seqnumpair_set1.length + seqnumpair_set2.length
  size_both = (seqnumpair_set1 & seqnumpair_set2).length
  perc_both = showcomment(2 * size_both,sum_size,"occur in set 1 and set 2")
  size_only_1 = (seqnumpair_set1 - seqnumpair_set2).length
  if size_only_1 > 0
    perc_only_1 = showcomment(size_only_1,sum_size,
                                   "occur in set 1 but not set 2")
  else
    perc_only_1 = 0
  end
  size_only_2 = (seqnumpair_set2 - seqnumpair_set1).length
  if size_only_2 > 0
    perc_only_2 = showcomment(size_only_2,sum_size,
                                  "occur in set 2 but not set 1")
  else
    perc_only_2 = 0
  end
  return [sum_size,perc_both,perc_only_1,perc_only_2]
end

def fill_other_hash(h1,h2)
  h1.each_pair do |k,v1|
    if h2.has_key?(k)
      v2 = h2[k]
      if [v1.seq1,v1.seq2] != [v2.seq1,v2.seq2]
        STDERR.puts "seq: #{taglist[0]}=[#{v1.seq1},#{v1.seq2}] != " +
                    "     [#{v2.seq1},#{v2.seq2}]=#{taglist[1]}"
        exit 1
      end
      h1[k].result1 = evaluate_itv(v1.start1,v1.len1,v2.start1,v2.len1)
      h1[k].result2 = evaluate_itv(v1.start2,v1.len2,v2.start2,v2.len2)
      h1[k].other = v2
    end
  end
end

def percent(a,b)
  return 100.0 * a.to_f/b.to_f
end

def cmpseedhashes(checkbetter,minidentity,taglist,h1,h2,silent = false)
  fill_other_hash(h1,h2)
  nobrother = 0
  h1.each_pair do |k,v1|
    if v1.identity >= minidentity and not h2.has_key?(k)
      if not silent
        puts "#{taglist[0]}: #{k}=>#{match_to_s(v1)}, #{taglist[1]}=[]"
      end
      nobrother += 1
    end
  end
  count_identical = 0
  count_contained = 0
  count_containing = 0
  ovperc_hash = Hash.new()
  miniddiff = 3
  commentcheckbetter = "matches in #{taglist[1]} >= #{miniddiff}%-points " +
                       "better than the corresponding matches in " +
                       "#{taglist[0]}:"
  h1.each_pair do |k,v|
    if v.result2.nil?
      next
    end
    if (v.result1.contained and v.result2.contained) or
       (v.result1.containing and v.result2.containing)
      count_identical += 1
    elsif v.result1.contained or v.result2.contained
      count_contained += 1
    elsif v.result1.containing or v.result2.containing
      count_containing += 1
    elsif checkbetter
      if v.other.identity >= minidentity
        if (v.other.len1 > v.len1 or v.other.len2 > v.len2) and
           (v.other.identity > v.identity) and 
           (v.other.identity - v.identity >= miniddiff.to_f)
          if not commentcheckbetter.nil?
             puts commentcheckbetter
             commentcheckbetter = nil
          end
          puts "#{taglist[0]}=#{match_to_s(v)}"
          puts "#{taglist[1]}=#{match_to_s(v.other)}"
          addpercentage(ovperc_hash,v.result1.ovperc)
          addpercentage(ovperc_hash,v.result2.ovperc)
        end
      end
    end
  end
  h1perc = percent(count_identical,h1.length)
  h2perc = percent(count_identical,h2.length)
  puts "matches in #{taglist[0]}: #{h1.length}"
  puts "matches in #{taglist[1]}: #{h2.length}"
  printf("identical=%d: %.2f of %s/%.2f of %s\n",count_identical,
             h1perc,taglist[0],h2perc,taglist[1])
  puts "contained=#{count_contained}"
  puts "containing=#{count_containing}"
  ovperc_hash.sort.each do |k|
    puts "#{k[0]}%\t#{k[1]}"
  end
  puts "matches in #{taglist[0]} with no corresponding matches in " +
       "#{taglist[1]}=#{nobrother}"
end

def cmpextendedlength(minidentity,h1,h2)
  fill_other_hash(h1,h2)
  lendiff_dist = Hash.new()
  h1.each_pair do |k,v1|
    if v1.identity < minidentity
      next
    end
    h1len = (v1.len1 + v1.len2)/2
    if not h2.has_key?(k)
      h2len = 0
    else
      h2len = (h2[k].len1 + h2[k].len2)/2
    end
    lendiff = (100.0 * (h1len - h2len).to_f/[h1len,h2len].max.to_f).to_i
    if not lendiff_dist.has_key?(lendiff)
      lendiff_dist[lendiff] = 1
    else
      lendiff_dist[lendiff] += 1
    end
  end
  h2.each_pair do |k,v2|
    if v2.identity < minidentity
      next
    end
    h2len = (v2.len1 + v2.len2)/2
    if not h1.has_key?(k)
      lendiff = -100
      if not lendiff_dist.has_key?(lendiff)
        lendiff_dist[lendiff] = 1
      else
        lendiff_dist[lendiff] += 1
      end
    end
  end
  lendiff_dist.sort.each do |diff,count|
    puts "diff #{diff}\t#{count}"
  end
end

def preprocess_index(inputfile)
  suffixeratorcall = "env -i bin/gt suffixerator -suftabuint " +
                     "-db #{inputfile} " +
                     "-dna -suf -tis -lcp -md5 no -des no -sds no "
  indexname = File.basename(inputfile)
  if not File.exist?("#{indexname}.prj") or
    File.stat("#{indexname}.prj").mtime < File.stat(inputfile).mtime
    puts "# #{suffixeratorcall}"
    if not system(suffixeratorcall)
      STDERR.puts "FAILURE: #{suffixeratorcall}"
      exit 1
    end
  end
  return indexname
end

def readmatchesfromfile(resultmatchfile)
  thishash = Hash.new()
  key = 0
  File.open(resultmatchfile,"r").each_line do |line|
    if line.match(/^\d/)
      thishash[key] = match_new_without_result(line.split(/\s/))
      key += 1
    elsif not line.match(/^\#/)
      STDERR.print "#{$0}: #{resultfile}: illegal line #{line}"
      exit 1
    end
  end
  return thishash
end
