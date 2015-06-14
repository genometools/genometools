#!/usr/bin/env ruby

require "set"
require 'ostruct'
require 'optparse'

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

def makeseedhash(indexname,seedlength,errperc,maxalilendiffopt,
                 extend_opt)
  seedhash = Hash.new()
  key = nil
  if seedlength == 0
    seedlength = 3 * prefixlength_get("#{indexname}.prj")
  end
  repfindcall = "env -i bin/gt repfind -scan -v -seedlength #{seedlength} " +
                "-#{extend_opt} -ii #{indexname} -err #{errperc}" +
                " #{maxalilendiffopt}"
  puts "\# #{repfindcall}"
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
      seedhash[key] = Match.new(a[0].to_i,a[1].to_i,a[2].to_i,
                                a[4].to_i,a[5].to_i,a[6].to_i,
                                a[7].to_i,a[8].to_i,a[9].to_f,
                                nil,nil,nil)
      key = nil
    end
  end
  if "#{$?}" != "" and not "#{$?}".match(/exit 0$/)
    STDERR.puts "FAILURE: #{repfindcall}: \"#{$?}\""
    exit 1
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
    printf("# %d (%.0f%%) of %d sequence pairs %s\n",size,perc,sum_size,comment)
    return perc.to_i
  else 
    return 0
  end
end

def calcdifference(seqnumpair_set1,seqnumpair_set2)
  sum_size = seqnumpair_set1.length + seqnumpair_set2.length
  size_both = (seqnumpair_set1 & seqnumpair_set2).length
  perc_both = showcomment(size_both,sum_size,"occur in greedy and xdrop")
  size_only_greedy = (seqnumpair_set1 - seqnumpair_set2).length
  perc_only_greedy = showcomment(size_only_greedy,sum_size,
                                 "occur in greedy but not xdrop")
  size_only_xdrop = (seqnumpair_set2 - seqnumpair_set1).length
  perc_only_xdrop = showcomment(size_only_xdrop,sum_size,
                                "occur in xdrop but not greedy")
  return [sum_size,perc_both,perc_only_greedy,perc_only_xdrop].join("\t")
end

def cmpseedhashes(checkbetter,minidentity,taglist,h1,h2)
  nobrother = 0
  h1.each_pair do |k,v1|
    if not h2.has_key?(k)
      if v1.identity >= minidentity
        puts "#{taglist[0]}: #{k}=>#{match_to_s(v1)}, #{taglist[1]}=[]"
        nobrother += 1
      end
    else
      v2 = h2[k]
      if [v1.seq1,v1.seq2] != [v2.seq1,v2.seq2]
        puts "seq: #{taglist[0]}=[#{v1.seq1},#{v1.seq2}] != " +
             "     [#{v2.seq1},#{v2.seq2}]=#{taglist[1]}"
        exit 1
      end
      h1[k].result1 = evaluate_itv(v1.start1,v1.len1,v2.start1,v2.len1)
      h1[k].result2 = evaluate_itv(v1.start2,v1.len2,v2.start2,v2.len2)
      h1[k].other = v2
    end
  end
  count_identical = 0
  count_contained = 0
  count_containing = 0
  ovperc_hash = Hash.new()
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
           (v.other.identity > v.identity)
          puts "#{taglist[0]}=#{match_to_s(v)}"
          puts "#{taglist[1]}=#{match_to_s(v.other)}"
          addpercentage(ovperc_hash,v.result1.ovperc)
          addpercentage(ovperc_hash,v.result2.ovperc)
        end
      end
    end
  end
  puts "identical=#{count_identical}"
  puts "contained=#{count_contained}"
  puts "containing=#{count_containing}"
  ovperc_hash.sort.each do |k|
    puts "#{k[0]}%\t#{k[1]}"
  end
  puts "#{taglist[0]}: nobrother=#{nobrother}"
end

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.inputfile = nil
  options.seedlength = 20
  options.errperc = 10
  options.maxalilendiff = 30
  options.silent = false
  opts = OptionParser.new
  opts.on("--inputfile STRING","specify input file") do |x|
    options.inputfile = x
  end
  opts.on("--seedlength NUM","specify seedlength") do |x|
    options.seedlength = x.to_i
  end
  opts.on("--errperc NUM","specify error percentage") do |x|
    options.errperc = x.to_i
  end
  opts.on("--maxalilendiff NUM",
          "specify maximum difference of aligned length") do |x|
    options.maxalilendiff = x.to_i
  end
  opts.on("--silent","do not show differences") do |x|
    options.silent = true
  end
  rest = opts.parse(argv)
  if not rest.empty?
    STDERR.puts "#{$0}: superfluous arguments: #{rest}"
    exit 1
  end
  if options.inputfile.nil?
    STDERR.puts "#{$0}: option -inputfile is mandatory"
    exit 1
  end
  return options
end

options = parseargs(ARGV)
  
suffixeratorcall = "env -i bin/gt suffixerator -suftabuint " +
                   "-db #{options.inputfile} " +
                   "-dna -suf -tis -lcp -md5 no -des no -sds no "
indexname = File.basename(options.inputfile)
if not File.exist?("#{indexname}.prj") or 
  File.stat("#{indexname}.prj").mtime < File.stat(options.inputfile).mtime
  if not system(suffixeratorcall)
    STDERR.puts "FAILURE: #{suffixeratorcall}"
    exit 1
  end
  puts "# #{suffixeratorcall}"
end
taglist = ["greedy","xdrop"]
seedhash1 = makeseedhash(indexname,options.seedlength,options.errperc,
                         "-maxalilendiff #{options.maxalilendiff}",
                         "extend#{taglist[0]}")
puts "# seedhash1: size = #{seedhash1.length}"
seedhash2 = makeseedhash(indexname,options.seedlength,options.errperc,"",
                         "extend#{taglist[1]}")
puts "# seedhash2: size = #{seedhash2.length}"

minidentity = 100 - options.errperc
seqnumpair_set1 = seedhash2seqnum_pairs(seedhash1)
seqnumpair_set2 = seedhash2seqnum_pairs(seedhash2)
result = calcdifference(seqnumpair_set1,seqnumpair_set2)
puts "#{options.seedlength}\t#{result}"
if not options.silent
  cmpseedhashes(true,minidentity,taglist,seedhash1,seedhash2)
  cmpseedhashes(false,minidentity,taglist.reverse,seedhash2,seedhash1)
end
