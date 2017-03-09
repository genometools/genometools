#!/usr/bin/env ruby

require "set"
require "optparse"

SEmatch = Struct.new("SEmatch",:dbseqnum,
                               :dbstart,
                               :dbend,
                               :queryseqnum,
                               :querystart,
                               :queryend,
                               :forward,
                               :score,
                               :origline)

def openfile(filename)
begin
  fp = File.new(filename,"r")
rescue => err
  STDERR.puts "#{$0}: cannot open #{filename}: #{err}"
end
  return fp
end

def coords_contained_in(start0,end0,start1,end1)
  if start1 <= start0 and end0 <= end1
    return true
  end
  return false
end

def match_is_identical(m0,m1)
  if m0.dbseqnum != m1.dbseqnum or m0.queryseqnum != m1.queryseqnum
    STDERR.puts "expect same sequence numbers"
    exit 1
  end
  if [m0.dbstart,m0.dbend,m0.querystart,m0.queryend] ==
     [m1.dbstart,m1.dbend,m1.querystart,m1.queryend]
    return true
  end
  return false
end

def match_proper_contained_in(m0,m1)
  if m0.dbseqnum != m1.dbseqnum or m0.queryseqnum != m1.queryseqnum
    STDERR.puts "expect same sequence numbers"
    exit 1
  end
  if coords_contained_in(m0.dbstart,m0.dbend,m1.dbstart,m1.dbend) and
     coords_contained_in(m0.querystart,m0.queryend,m1.querystart,m1.queryend) and
     not match_is_identical(m0,m1)
    return true
  end
  return false
end

def coords_overlap_size(start0,end0,start1,end1)
  if start0 > start1
    STDERR.puts "start0=#{start0} > #{start1} = start1 not expected"
    exit 1
  end
  if end0 < start1
    return 0
  elsif end0 < end1
    return end0 - start1 + 1
  else
    return end1 - start1 + 1
  end
end

def matchlength(m)
  return (m.dbend - m.dbstart + 1 + m.queryend - m.querystart + 1)/2
end

def matches_overlap(m0,m1)
  ovl = 0
  if m0.dbstart <= m1.dbstart
    ovl += coords_overlap_size(m0.dbstart,m0.dbend,m1.dbstart,m1.dbend)
  else
    ovl += coords_overlap_size(m1.dbstart,m1.dbend,m0.dbstart,m0.dbend)
  end
  if m0.querystart <= m1.querystart
    ovl += coords_overlap_size(m0.querystart,m0.queryend,m1.querystart,m1.queryend)
  else
    ovl += coords_overlap_size(m1.querystart,m1.queryend,m0.querystart,m0.queryend)
  end
  len = matchlength(m0) + matchlength(m1)
  return (ovl.to_f/len.to_f).round
end

def add_single_match(matchset,m,always_add = false)
  add_m = true
  if not always_add
    to_delete = Array.new()
    matchset.each do |previous|
      if match_proper_contained_in(m,previous)
        if previous.score > m.score
          add_m = false
        end
      elsif match_proper_contained_in(previous,m)
        if previous.score < m.score
          to_delete.push(previous)
        end
      elsif match_is_identical(previous,m) and previous.score > m.score
        add_m = false
      end
    end
    to_delete.each do |elem|
      matchset.delete(elem)
    end
  end
  if add_m
    matchset.add(m)
  end
  max_score = 0
  matchset.each do |elem|
    if max_score < elem.score
      max_score = elem.score
    end
  end
  return max_score
end

def convertmatchfile2hash(matchfile,always_add)
  match_hash = Hash.new()
  score_hash = Hash.new() {0}
  fp = openfile(matchfile)
  runtime = nil
  spacepeak = nil
  fp.each_line do |line|
    if line.match(/^#/)
      m = line.match(/^# TIME.*\s(\S+)$/)
      if m
        runtime = m[1].to_f
      elsif m = line.match(/^# combined space peak in megabytes:\s(\S+)$/)
        spacepeak = m[1].to_f
      end
    else
      a = line.split(/\s/)
      m = SEmatch.new(a[1].to_i,
                      a[2].to_i,
                      a[2].to_i + a[0].to_i - 1,
                      a[5].to_i,
                      a[6].to_i,
                      a[6].to_i + a[4].to_i - 1,
                      if a[3] == "F" then true else false end,
                      a[7].to_i,
                      line.chomp)
      key = [m.dbseqnum,m.queryseqnum]
      if not match_hash.has_key?(key)
        match_hash[key] = Set.new()
      end
      score_hash[key] = add_single_match(match_hash[key],m,always_add)
    end
  end
  return match_hash, score_hash, runtime, spacepeak
end

def merge_sets(mh0,mh1)
  seqnumpairkeys = (mh0.keys + mh1.keys).to_set
  seqnumpairkeys.each do |key|
    mh0_set = if mh0.has_key?(key) then mh0[key] else Set.new() end
    mh1_set = if mh1.has_key?(key) then mh1[key] else Set.new() end
    yield key, mh0_set, mh1_set
  end
end

def multi_merge_sets(mh_list)
  seqnumpairkeys = Set.new()
  mh_list.each do |mh|
    mh.keys.each do |key|
      seqnumpairkeys.add(key)
    end
  end
  seqnumpairkeys.each do |key|
    mh_sets = Array.new(mh_list.length)
    mh_list.each_with_index do |mh,idx|
      mh_sets[idx] = if mh.has_key?(key) then mh[key] else Set.new() end
    end
    yield key, mh_sets
  end
end

Counters = Struct.new("Counters",:different_seqpairs,
                                 :different_matches,
                                 :identical_seqpairs,
                                 :identical_matches,
                                 :mh0larger,
                                 :mh1larger,
                                 :samesize,
                                 :mh0_all_overlap,
                                 :mh0_not_all_overlap,
                                 :mh1_all_overlap,
                                 :mh1_not_all_overlap,
                                 :mh0_largerscore,
                                 :mh1_largerscore,
                                 :samescore)

def overlap_in_one_instance(mh0_set,mh1_set)
  count_overlaps = 0
  mh0_set.each do |m0|
    overlaps = false
    mh1_set.each do |m1|
      if coords_contained_in(m0.dbstart,m0.dbend,m1.dbstart,m1.dbend) or
         coords_contained_in(m0.querystart,m0.queryend,
                             m1.querystart,m1.queryend)
        overlaps = true
        break
      end
    end
    if overlaps
      count_overlaps += 1
    end
  end
  return count_overlaps
end

def set_show(key,s)
  if s.empty?
    puts "# #{key} is empty"
  else
    s.each do |elem|
      puts "# #{key}: #{elem.origline}"
    end
  end
end

def hash_difference(mh0,sh0,mh1,sh1)
  counters = Counters.new(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  multi_merge_sets([mh0,mh1]) do |key, mh_sets|
    mh0_set = mh_sets[0]
    mh1_set = mh_sets[1]
    if mh0_set == mh1_set
      counters.identical_matches += mh0_set.size
      counters.identical_seqpairs += 1
    else
      counters.different_seqpairs += 1
      common = mh0_set.intersection(mh1_set)
      counters.identical_matches += common.length
      mh0_set_alone = mh0_set.difference(common)
      mh1_set_alone = mh1_set.difference(common)
      counters.different_matches += mh0_set_alone.size + mh1_set_alone.size
      set_show("common ",common)
      set_show("mh0_set minus common",mh0_set_alone)
      set_show("mh1_set minus common",mh1_set_alone)
      puts ""
      if mh0_set.size > mh1_set.size
        counters.mh0larger += 1
      elsif mh0_set.size < mh1_set.size
        counters.mh1larger += 1
      else
        counters.samesize += 1
      end
      overlaps = overlap_in_one_instance(mh1_set,mh0_set)
      if overlaps == mh1_set.size
        counters.mh1_all_overlap += 1
      else
        counters.mh1_not_all_overlap += 1
      end
      overlaps = overlap_in_one_instance(mh0_set,mh1_set)
      if overlaps == mh0_set.size
        counters.mh0_all_overlap += 1
      else
        counters.mh0_not_all_overlap += 1
      end
      if sh0[key] > sh1[key]
        if sh0[key] > sh1[key] * 11/10
          counters.mh0_largerscore += 1
        else
          counters.samescore += 1
        end
      elsif sh0[key] < sh1[key]
        if sh0[key] * 11/10 < sh1[key]
          counters.mh1_largerscore += 1
        else
          counters.samescore += 1
        end
      else
        counters.samescore += 1
      end
    end
  end
  return counters
end

def identical_upto_query(s,t)
  if s.size != t.size
    STDERR.puts "#{$0} are of different sizes #{s.size} and #{t.size}"
    exit 1
  end
  s.zip(t).each do |ms,mt|
    vals = [ms.dbseqnum,ms.dbstart,ms.dbend,ms.queryend-ms.querystart,ms.score,
            ms.forward]
    valt = [mt.dbseqnum,mt.dbstart,mt.dbend,mt.queryend-mt.querystart,mt.score,
            if mt.forward then false else true end]
    if vals != valt
      STDERR.puts "#{$0}: vals = #{vals} != #{valt}"
      exit 1
    end
  end
end

MCoptions = Struct.new("MCoptions",:pairwise,:ignore_query,:matchfiles)

def parseargs(argv)
  options = MCoptions.new(false,false,nil)
  opts = OptionParser.new
  opts.banner = "#{$0} [options] <inputfile>"
  opts.on("-p","--pairwise","perform pairwise comparison") do |x|
    options.pairwise = true
  end
  opts.on("-i","--ignore-query-positions","comparison, ignoring sequence " +
                                          "number and start position " +
                                          "on query") do |x|
    options.ignore_query = true
  end
  rest = opts.parse(argv)
  if rest.length < 2
    STDERR.puts "Usage: #{$0}: missing input files"
    exit 1
  else
    options.matchfiles = rest
  end
  if options.pairwise and options.ignore_query
    STDERR.puts "#{$0}: options -p and -i exclude each other"
    exit 1
  end
  return options
end

options = parseargs(ARGV)

numfiles = options.matchfiles.length
match_hash_tab = Array.new(numfiles)
score_hash_tab = Array.new(numfiles)
runtime_tab = Array.new(numfiles)
spacepeak_tab = Array.new(numfiles)
always_add = if options.ignore_query then true else false end
0.upto(numfiles-1).each do |idx|
  match_hash_tab[idx], score_hash_tab[idx], runtime_tab[idx],
  spacepeak_tab[idx] = convertmatchfile2hash(options.matchfiles[idx],always_add)
end

if options.pairwise
  0.upto(numfiles-2).each do |idx0|
    argv0 = options.matchfiles[idx0]
    (idx0+1).upto(numfiles-1).each do |idx1|
      argv1 = options.matchfiles[idx1]
      puts "# #{argv0} vs #{argv1}"
      counters = hash_difference(match_hash_tab[idx0],score_hash_tab[idx0],
				 match_hash_tab[idx1],score_hash_tab[idx1])
      all_seqpairs = counters.identical_seqpairs + counters.different_seqpairs
      printf("sequence pairs with identical matches %d (%5.2f%%)\n",
	    counters.identical_seqpairs,
	    100.0 * counters.identical_seqpairs.to_f/all_seqpairs.to_f)
      printf("sequence pairs with different matches %d (%5.2f%%)\n",
	      counters.different_seqpairs,
	      100.0 * counters.different_seqpairs.to_f/all_seqpairs.to_f)
      all_matches = counters.identical_matches + counters.different_matches
      printf("identical matches %d (%5.2f%%)\n",
	      counters.identical_matches,
	      100.0 * counters.identical_matches.to_f/all_matches.to_f)
      printf("different matches %d (%5.2f%%)\n",
	      counters.different_matches,
	      100.0 * counters.different_matches.to_f/all_matches.to_f)

      puts "sequence pairs for which #{argv0} > #{argv1}: #{counters.mh0larger}"
      puts "sequence pairs for which #{argv1} > #{argv0}: #{counters.mh1larger}"
      puts "sequence pairs for which #{argv1} = #{argv0}: #{counters.samesize}"

      puts "sequence pairs with different match sets for which all of #{argv0} overlaps with #{argv1}: #{counters.mh0_all_overlap}"
      puts "sequence pairs with different match sets for which some of #{argv0} do not overlap with #{argv1}: #{counters.mh0_not_all_overlap}"

      puts "sequence pairs with different match sets for which all of #{argv1} overlaps with #{argv0}: #{counters.mh1_all_overlap}"
      puts "sequence pairs with different match sets for which some of #{argv1} do not overlaps with #{argv0}: #{counters.mh1_not_all_overlap}"

      puts "sequence pairs with different match sets for which #{argv1} have a larger score than #{argv0}: #{counters.mh1_largerscore}"
      puts "sequence pairs with different match sets for which #{argv0} have a larger score than #{argv1}: #{counters.mh0_largerscore}"
      puts "sequence pairs for which #{argv0} have same score as #{argv1}: #{counters.samescore}"
    end
  end
elsif options.ignore_query
  multi_merge_sets(match_hash_tab) do |key,mh_sets|
    mh_sets.each_with_index do |s,sidx|
      mh_sets.each_with_index do |t,tidx|
        if sidx < tidx
          identical_upto_query(s,t)
        end
      end
    end
  end
else
  contained = Array.new(numfiles) {0}
  step = 10
  bestscore = Array.new(numfiles) {Array.new(100/step+1) {0}}
  totalsetsize = 0
  numsets = 0
  mh_set_combined = Set.new()
  multi_merge_sets(match_hash_tab) do |key,mh_sets|
    max_score = 0
    mh_sets.each do |mh_set|
      mh_set.each do |m|
        max_score = add_single_match(mh_set_combined,m)
      end
    end
    mh_sets.each_with_index do |mh_set,midx|
      inse = mh_set_combined.intersection(mh_set)
      contained[midx] += inse.size
      if score_hash_tab[midx][key] > max_score
        STDERR.puts "midx #{midx}: score_hash[#{key}]=" +
                    "#{score_hash_tab[midx][key]} > #{max_score}=max_score"
        exit(1)
      end
      if score_hash_tab[midx][key] > 0
        diff = (100 * (max_score - score_hash_tab[midx][key]).to_f/max_score.to_f).round
        bestscore[midx][diff/step] += 1
      else
        bestscore[midx][100/step] += 1
      end
    end
    numsets += 1
    totalsetsize += mh_set_combined.size
    mh_set_combined.clear()
  end
  result = Array.new()
  length_of_longest = 0
  0.upto(numfiles-1).each do |idx|
    result.push([bestscore[idx][0],idx])
    if length_of_longest < options.matchfiles[idx].length
      length_of_longest = options.matchfiles[idx].length
    end
  end
  column_header = Range.new(step,100).step(step).to_a.map {|x| if x < 100 then " lt#{x}" else "lt#{x}" end}
  printf("%*s\t%5s\t%s\t  nil\t time\tspace\n",length_of_longest,"options","cont",column_header.join("\t"))
  result.sort.each do |cont,idx|
    printf("%*s\t%5.2f",length_of_longest,options.matchfiles[idx],
                        100.0 * cont.to_f/totalsetsize.to_f)
    sumpercentage = 0.0
    0.upto(100/step-1).each do |diff|
      bestpercentage = 100.0 * bestscore[idx][diff].to_f/numsets.to_f
      sumpercentage += bestpercentage
      printf("\t%5.2f",sumpercentage)
    end
    printf("\t%5.2f\t%.2f\t%.2f",
                100.0 * bestscore[idx][100/step].to_f/numsets.to_f,
                runtime_tab[idx],spacepeak_tab[idx])
    puts ""
  end
end
