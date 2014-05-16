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
  options.amalgamate = false
  options.ccache = true
  options.check = true
  options.clang = true
  options.ddd = false
  options.gprof = false
  options.m64 = false
  options.prof = false
  options.speed = false
  options.makethreads = false
  options.threds = false
  opts = OptionParser.new
  opts.on("-a", "--amalgamation",
          "Combine code into one file for optimisation") do |x|
    options.amalgamate = true
  end
  opts.on("-c", "--[no-]clang", "use [no] clang, defaults to using it") do |x|
    options.clang = x
  end
  opts.on("-d","--ddd","compile for debugging (no opt)") do |x|
    if options.speed
      usage(opts,"-d and -s exclusive")
    end
    options.ddd = true
  end
  opts.on("-f", "--force", "do not check if make is the same as last") do |x|
    options.check = false
  end
  opts.on("-g","--gprof","compile for google-proftool") do |x|
    options.gprof = true
    if options.prof
      usage(opts,"-p and -g are exclusive")
    end
  end
  opts.on("-j [n]", Integer, "number of processes for make") do |x|
    options.makethreads = true
    options.j = x || ""
  end
  opts.on("-m","--m64","compile 64 bit binary") do |x|
    options.m64 = true
  end
  opts.on("-p","--prof","compile for profiling") do |x|
    options.prof = true
    if options.gprof
      usage(opts,"-p and -g are exclusive")
    end
  end
  opts.on("-s","--speed","optimize for speed, implies -a") do |x|
    if options.ddd
      usage(opts,"-d and -s exclusive")
    end
    options.speed = true
  end
  opts.on("-t", "--threads", "compile with threading enabled") do |x|
    options.check = false
  end
  opts.on("--[no-]ccache", "use [no] ccache, defaults to using it") do |x|
    options.ccache = x
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"unnecessary arguments: #{rest}")
    exit 1
  end
  puts options
  return options
end

def makecompilerflags(fp,options)
  fp.print "all:\n\t\${MAKE} cairo=no"
  fp.print " CFLAGS+=-fstrict-aliasing"
  if options.ddd
    fp.print " opt=no"
  end
  if options.speed
    fp.print " assert=no amalgamation=yes"
  elsif options.amalgamate
    fp.print " amalgamation=yes"
  end
  if options.m64
    fp.print " 64bit=yes"
  end
  if options.prof
    fp.print " prof=yes"
    options.clang = false
  end
  if options.ccache
    if options.clang
      fp.print " CC='ccache clang'"
    else
      fp.print " CC='ccache gcc'"
    end
  else
    if options.clang
      fp.print " CC='clang'"
    else
      fp.print " CC='gcc'"
    end
  end
  if options.gprof
    fp.print " LIBS='-Wl,--no-as-needed -lprofiler -Wl,--as-needed'"
  end
  if options.threads
    fp.print " threads=yes"
  end
  if options.makethreads
    fp.print " -j#{options.j}"
  end
  fp.puts
end

if File.exists?('LocalMakefile')
  FileUtils.mv('LocalMakefile','LocalMakefile.previous')
end

options = parseargs(ARGV)

File.open('LocalMakefile',"w") do |fp|
  makecompilerflags(fp,options)
end

if File.exists?('LocalMakefile.previous') and
    not FileUtils.compare_file('LocalMakefile','LocalMakefile.previous') and
    options.check
  STDERR.puts "Current and previous LocalMakefile files differ: first " +
              "remove them"
  exit 1
end

system("make -f LocalMakefile")
