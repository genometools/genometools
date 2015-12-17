#!/usr/bin/env ruby

require 'optparse'
require 'ostruct'

Sortinfo = Struct.new("Sortinfo",:len,:seqnum,:startfwdstrand)

def line2sortinfo(line)
  arr = line.split(/\s/)
  return Sortinfo.new(arr[4].to_i,arr[5].to_i,arr[6].to_i)
end

def compare_matches(line0,line1)
  sortinfo0 = line2sortinfo(line0)
  sortinfo1 = line2sortinfo(line1)
  if sortinfo0.seqnum < sortinfo1.seqnum or
       (sortinfo0.seqnum == sortinfo1.seqnum and
        sortinfo0.startfwdstrand + sortinfo0.len <=
        sortinfo1.startfwdstrand + sortinfo1.len)
    return -1
  end
  return 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.inputfile = nil
  opts = OptionParser.new()
  opts.on("-c","--check","check order") do |x|
    options.check = true
  end
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts "Usage: #{$0} [options] [inputfile]"
    puts "read matches from file or STDIN and sort them in ascending order of\nquery numbers. Matches on the same query are sorted in ascending order of the end position in query"
    exit 0
   end
  rest = opts.parse(argv)
  if rest.length == 1
    options.inputfile = rest[0]
  elsif rest.length > 1
    STDERR.puts options.banner
    exit 1
  end
  return options
end

def matches_checkorder(fpin)
  previousline = nil
  fpin.each_line do |line|
    if line.match(/^[0-9]/)
      if not previousline.nil? and compare_matches(previousline,line) > 0
        STDERR.puts "#{previousline.chomp} >\n#{line}"
        exit 1
      end
      previousline = line
    end
  end
end

def matches_orderedoutput(fpin)
  header = Array.new()
  headercomplete = false
  footer = Array.new()
  matchlines = Array.new()
  fpin.each_line do |line|
    if line.match(/^#/)
      if not line.match(/^# seed/)
        if headercomplete
          footer.push(line)
        else
          header.push(line)
        end
      end
    elsif line.match(/^[0-9]/)
      headercomplete = true
      matchlines.push(line)
    end
  end
  matchlines.sort! {|a,b| compare_matches(a,b)}
  puts header
  puts matchlines
  puts footer
end


options = parseargs(ARGV)
if options.inputfile.nil?
  fpin = STDIN
else
  begin
    fpin = File.new(options.inputfile,"r")
  rescue => err
    STDERR.puts "#{$0}: cannot open #{options.inputfile}: #{err}"
    exit 1
  end
end
if options.check
  matches_checkorder(fpin)
else
  matches_orderedoutput(fpin)
end
if not options.inputfile.nil?
  fpin.close_read
end
