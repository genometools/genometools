#!/usr/bin/env ruby

require_relative "turnwheel.rb"

elements = ["-tis","-suf","-bwt","-lcp"]
size = elements.size
turnwheels(Array.new(size) {2}) do |wheel|
  0.upto(size-1).each do |i|
    if wheel[i] == 0
      print " #{elements[i]}"
    end
  end
  puts ""
end
