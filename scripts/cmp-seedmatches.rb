#!/usr/bin/env ruby

require "scripts/evalseedhash.rb"

def filename2keys(filename)
  File.open(filename).each_line do |line|
    if m = line.match(/^\# TIME [A-Za-z]+-([a-z]+)*([\d\-]+)/)
      minidentity = (line.split(/-/)[4]).to_i
      return m[1], m[2], minidentity
    end
  end
  return nil
end

seedhash_tab = Hash.new()
filelist = Array.new()
ARGV.each do |matchfile|
  STDERR.print "read #{matchfile} "
  seedhash_tab[matchfile] = inputseedhash(matchfile)
  STDERR.puts "contains #{seedhash_tab[matchfile].length} matches"
  filelist.push(matchfile)
end

0.upto(filelist.length-1) do |idx1|
  filename1 = filelist[idx1]
  tag1, spec1, minidentity1 = filename2keys(filename1)
  0.upto(filelist.length-1) do |idx2|
    filename2 = filelist[idx2]
    tag2, spec2, minidentity2 = filename2keys(filename2)
    if idx1 < idx2
      puts "minidentity=#{minidentity1}%: #{tag1}#{spec1} vs #{tag2}#{spec2}"
      cmpseedhashes(true,100.0 - minidentity1,[tag1,tag2],
                    seedhash_tab[filename1],seedhash_tab[filename2],true)
      cmpseedhashes(false,100.0 - minidentity1,[tag2,tag1],
                    seedhash_tab[filename2],seedhash_tab[filename1],true)
    end
  end
end
