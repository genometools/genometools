#!/usr/bin/env ruby
# compare two files of tabular output by hmmsearch and calculate FP NP etc.
# only needs to get the first element of each row, as whole sequences are
# reported

require 'set'

testset = Hash.new{|hash,key| hash[key] = Set.new}
fpset = Hash.new{|hash,key| hash[key] = Set.new}
tpset = Hash.new{|hash,key| hash[key] = Set.new}

File.open(ARGV[0]) do |file|
  file.each_line do |line|
    next if line[0].match /^\s*#/
    if match = line.match /^(\S+)\s+(\S+)/
      testset[match[2]].add match[1]
    else
      STDERR.puts: "couldn't parse: #{line}"
      exit 1
    end
  end
end

File.open(ARGV[1]) do |file|
  file.each_line do |line|
    next if line[0].match /^\s*#/
    if match = line.match /^(\S+)\s+(\S+)/
      if tpset[match[1]].member? match[0]
        tpset[match[1]].add match[0]
      else
        fpset[match[1]].add match[0]
      end
    else
      STDERR.puts: "couldn't parse: #{line}"
      exit 1
    end
  end
end

total = 0
tp = 0
fp = 0
fn = 0

testset.each_key do |key|
  total += testset[key].size
end

tpset.each_key do |key|
  tp += tpset[key].size
end

fpset.each_key do |key|
  fp += fpset[key].size
end

puts "# testset\t#{total}"
puts "# tp:\t#{tp}"
puts "# fp:\t#{fp}"

fpset.each_key do |key|
  puts "# #{key}"
  fpset[key].each do |hit|
    puts hit
  end
end

testset.subtract(tpset)
puts "# fn:\t#{testset.size}"

testset.each do |hit|
  puts hit
end
