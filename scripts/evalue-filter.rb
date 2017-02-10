#!/usr/bin/env ruby

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <evalue-filter> <file with subset>"
  exit 1
end

def openthefile(filename)
  begin
    fp = File.new(filename,"r")
  rescue => err
    STDERR.puts "#{$0}: cannot open #{filename}: #{err}"
    exit 1
  end
end

evalue_filter = ARGV[0].to_f
fpstrong = openthefile(ARGV[1])

fpstrong.each_line do |line|
  a = line.split(/\s/)
  evalue = a[a.length - 1].to_f
  if evalue > evalue_filter
    STDERR.puts "#{$0}: unexpected line: #{line.chomp}"
    exit 1
  end
end
