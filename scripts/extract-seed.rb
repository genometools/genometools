#!/usr/bin/env ruby

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <inputfile>"
  exit 1
end

File.new(ARGV[0],"r").each_line do |line|
  if not line.match(/^#/)
    a = line.split(/\s/)
    if a.length < 13
      STDERR.puts "line #{line.chomp} contains #{a.length} < 13 columns"
      exit 1
    end
    puts "# seed: #{a[11]} #{a[12]} #{a[10]}"
  end
  print line
end
