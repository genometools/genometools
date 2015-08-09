#!/usr/bin/env ruby

require "set"
require "scripts/evalseedhash.rb"

def percdiff(a)
  return (a[5].to_i - a[6].to_i).abs
end

def scoreresult(a,maxdiff)
  total = a[3].to_i
  common = a[4].to_i
  return (common * common)/2 + total + 5 * (maxdiff - percdiff(a))
end

def maxdiff_get(matches)
  maxdiff = nil
  matches.each do |m|
    if maxdiff.nil? or maxdiff < percdiff(m)
      maxdiff = percdiff(m)
    end
  end
  return maxdiff
end

if ARGV.length != 4
  STDERR.puts "Usage: #{$0} <resultfile> <indexname> <seedlength> <minlength>"
  exit 1
end
resultfile = ARGV[0]
indexname = ARGV[1]
seedlength = ARGV[2].to_i
minlength = ARGV[3].to_i

len = nil
maxerrperc = 15
best = Array.new(maxerrperc+1) {nil}
bestsim = Array.new(maxerrperc+1) {Set.new()}
File.new(resultfile,"r").each_line do |line|
  a = line.chomp.split(/\t/)
  if not len.nil? and len != a.length
    STDERR.puts "lines of different length"
    exit 1
  else
    len = a.length
  end
  minidentity = a[0].to_i
  errperc = 100 - minidentity
  common = a[4].to_i
  if best[errperc].nil? or best[errperc][4].to_i < common
    best[errperc] = a
  end
end
File.new(ARGV[0],"r").each_line do |line|
  a = line.chomp.split(/\t/)
  minidentity = a[0].to_i
  errperc = 100 - minidentity
  common = a[4].to_i
  tolerance = (best[errperc][4].to_f/12.0).to_i
  if common + tolerance >= best[errperc][4].to_i
    bestsim[errperc].add(a)
  end
end

1.upto(maxerrperc).each do |errperc|
  maxdiff = maxdiff_get(bestsim[errperc])
  puts "------errperc=#{errperc}------maxdiff=#{maxdiff}"
  bestsim[errperc].sort {|x,y| scoreresult(y,maxdiff) <=> scoreresult(x,maxdiff)}.each do |a|
    puts (a.join("\t") + "\t" + scoreresult(a,maxdiff).to_s)
  end
end

puts "====== overall ====="
paramchoices = Array.new()
choosemin = 1
1.upto(maxerrperc).each do |errperc|
  maxdiff = maxdiff_get(bestsim[errperc])
  achoose = bestsim[errperc].select {|a| a[2].to_i >= choosemin}
  if achoose.nil?
    STDERR.puts "no entry with xdropbelow >= #{choosemin} in #{bestsim[errperc]}"
    exit
  end
  a = achoose.max_by {|x| scoreresult(x,maxdiff)}
  puts (a.join("\t") + "\t" + scoreresult(a,maxdiff).to_s)
  choosemin = a[2].to_i
  paramchoices[errperc] = [a[1].to_i,a[2].to_i]
end


emptyenv = false
1.upto(maxerrperc).each do |errperc|
  puts "#{paramchoices[errperc]}"
  commonoptions = ["-err #{errperc}","-l #{minlength}","-silent"]
  permathistory = paramchoices[errperc][0]
  xdropbelow = paramchoices[errperc][1]
  optionlist = commonoptions + ["-percmathistory #{permathistory}"]
  repfindcall = createrepfindcall(indexname,seedlength,"extendgreedy",
                                  optionlist,emptyenv,false)
  puts "#{repfindcall}"
  optionlist = commonoptions + ["-xdropbelow #{xdropbelow}"]
  repfindcall = createrepfindcall(indexname,seedlength,"extendxdrop",optionlist,
                                  emptyenv,false)
  puts "#{repfindcall}"
end
