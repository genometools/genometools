#!/usr/bin/env ruby
#
# Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

require 'ostruct'
require 'optparse'
require 'fileutils'

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.m64 = false
  options.speed = false
  options.prof = false
  opts = OptionParser.new
  opts.on("--m64","compile 64 bit binary") do |x|
    options.m64 = true
  end
  opts.on("--speed","optimize for speed") do |x|
    options.speed = true
  end
  opts.on("--prof","compile for profiling") do |x|
    options.prof = true
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"unnecessary arguments: #{rest}")
    exit 1
  end
  return options
end

def makecompilerflags(fp,options)
  fp.print "all:\n\t\${MAKE} curses=no cairo=no opt=yes"
  # fp.print " CFLAGS+=-DINLINEDSequentialsuffixarrayreader"
  if options.speed
    fp.print " assert=no amalgamation=yes"
  end
  if options.m64
    fp.print " 64bit=yes"
  end
  if options.prof
    fp.print " prof=yes"
  end
  fp.puts " CC='ccache gcc'"
end

if File.exists?('LocalMakefile')
  FileUtils.mv('LocalMakefile','LocalMakefile.previous')
end

options = parseargs(ARGV)

File.open('LocalMakefile',"w") do |fp|
  makecompilerflags(fp,options)
end

if File.exists?('LocalMakefile.previous') and
     not FileUtils.compare_file('LocalMakefile','LocalMakefile.previous')
  STDERR.puts "Current and previous LocalMakefile files differ: first " +
              "remove them"
  exit 1
end

system("make -f LocalMakefile")
