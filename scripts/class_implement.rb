#!/usr/bin/env ruby
#
# Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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
# # #
# creates an implementation skeleton of an abstract class
# see src/core/example_rep.h
#     src/core/example.[ch]
# for abstract class example and src/core/example_a.[ch]
#                                src/core/example_b.[ch]
# for implementations thereof.

require 'optparse'
require 'ostruct'
require 'pp'

require "erb"

$:.unshift File.join(File.dirname(__FILE__), ".")
require "codegen_module.rb"

module ClassGenOpts
  def self.parse(args)
    options = OpenStruct.new
    options.force = false

    opts = OptionParser.new()

    opts.program_name = File.basename($0)
    opts.version = 0 , 1
    opts.release = 1

    opts.banner = "USAGE: #{opts.program_name} [options] dirname classname basename"
    opts.separator ""
    opts.separator "dirname + classname gives the path to the abstract class,"
    opts.separator "basename will be added to classname for files of"
    opts.separator "implemetation"
    opts.separator ""
    opts.separator "Specific options:"

    opts.on("-f", "--force",
            "overwrite files if they already exist.") do |force|
      options.force = force
    end
    opts.separator ""
    opts.separator "Common options"
    opts.on("-h", "--help",
            "print this help end exit") do
      puts opts
      exit
    end
    opts.on("--version",
            "print version and exit") do
      puts opts.ver
      exit
    end

    opts.parse!(args)
    if (args.length != 3)
      CodeGen.perror "wrong number of arguments\n" + opts
    end
    options.dirname = File.join "src", args[0]
    options.classbase = File.join options.dirname, args[1]
    options.basename = options.classbase + '_' + args[2]

    if (File.exist?(options.basename + '.h') or
        File.exist?(options.basename + '.c')) and
        not options.force
      CodeGen.perror "Either or both file(s) #{options.basename}.[ch]\n" +
        "already exist. Use --force to overwrite.\n" +
        "(#{opts.program_name} -h for details)"
    end

    return options
  end
end

def create_header(options)
  h_file = File.new(options.basename + '.h', 'w')

  h_file << CodeGen.create_license

  classbase = options.classbase[4..-1]
  guard_macro = File.basename(options.basename).upcase + '_H'
  cvar= File.basename(options.classbase)
  icvar= File.basename(options.basename)
  fkt_pref = "gt_" + cvar
  ifkt_pref = "gt_" + icvar
  classN = CodeGen.filename2classname(options.classbase)
  iclassN = CodeGen.filename2classname(options.basename)

  content = ERB.new($implement_header).result(binding)

  h_file << ERB.new($hwrapper).result(binding)
  h_file.close
end

def create_src(options)
  src_file = File.new(options.basename + '.c', 'w')

  src_file << CodeGen.create_license

  classbase = options.classbase[/[^\/]+\/(.+)$/, 1]
  basename = options.basename[/[^\/]+\/(.+)$/, 1]
  classN = CodeGen.filename2classname(options.classbase)
  iclassN = CodeGen.filename2classname(options.basename)
  cvar= File.basename(options.classbase)
  icvar= File.basename(options.basename)
  fkt_pref = "gt_" + cvar
  ifkt_pref = "gt_" + icvar
  repfile = File.readlines(options.classbase + "_rep.h")
  functions = CodeGen.extract_functions_from_rep(classN,repfile)
  funcnames = []
  cvar= File.basename(options.classbase)
  icvar= File.basename(options.basename)
  src_file << ERB.new($implement_src).result(binding)
  src_file.close
end

begin
  options = ClassGenOpts.parse(ARGV);
rescue OptionParser::InvalidOption => e
  CodeGen.perror e.message
end

create_header(options)
create_src(options)
