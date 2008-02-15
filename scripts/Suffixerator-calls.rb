#!/usr/bin/env ruby

require 'Suffixerator-table'

numofargs = ARGV.length

if numofargs != 2
  STDERR.puts "Usage: $0 <program> <dbdna>"
  exit 1
end

program = ARGV[0]
dbdna = ARGV[1]

def makesystemcall(argstring)
  if system(argstring)
    STDERR.puts "# success #{argstring}"
  else
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1
  end
end

arglisttable = makesuffixeratorarglisttable(dbdna)

arglisttable.each do |argstring|
  makesystemcall(program + " " + argstring)
end

STDERR.puts "# number of calls was #{arglisttable.length}"
