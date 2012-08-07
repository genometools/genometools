#!/usr/bin/env ruby

ARGV.each do |filename|
  File.open(filename).each_line do |line|
    if line.match(/^>/)
      puts line
    else
      puts line.gsub(/[^ACGTacgt]/,"")
    end
  end
end
