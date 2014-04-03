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
# create a new abstract class, or uses an abstract class to start an
# implementation of it.
# see src/core/example_rep.h
#     src/core/example.[ch]
# for abstract class example and src/core/example_a.[ch]
#                                src/core/example_b.[ch]
# for implementations thereof.

require 'optparse'
require 'ostruct'

require "erb"

$:.unshift File.join(File.dirname(__FILE__), ".")
require "codegen_module.rb"

$CONFIG = ENV["GT_CODEGEN_CONFIG"] || "~/.gitconfig"

module ClassGenOpts
  def self.parse(args)
    options = OpenStruct.new
    options.editor = "editor"
    options.create = true
    options.force = false
    options.refcount = true

    opts = OptionParser.new()

    opts.program_name = File.basename($0)
    opts.version = 0 , 1
    opts.release = 1

    opts.banner = "USAGE: #{opts.program_name} [options] dirname basename"
    opts.separator ""
    opts.separator "editor will open during file-creation, after editing, safe"
    opts.separator "the file, more files will be created based on that file"
    opts.separator ""
    opts.separator "Specific options:"

    opts.on("--editor EDITOR",
            "Your editor of choice to edit source code, ",
            "defaults to 'editor' which is only set on ",
            "Debion style systems"
           ) do |editor|
      options.editor = editor
    end
    opts.on("-f", "--force",
            "overwrite files if they already exist. _rep.h ",
            "file can be reused with --no-create.") do |force|
      options.force = force
    end
    opts.on("--[no-]create",
            "defaults to --create, will overwrite _rep.h ",
            "file if it exists (with --force). Otherwise ",
            "it will try to reopen it for editing. ",
            "(usufull if something went wrong in the first ",
            "step)"
           ) do |create|
      options.create = create
    end
    opts.on("--[no-]refcount",
            "add reference counting to the class interface ",
            "defaults to true.") do |refcount|
      options.refcount = refcount
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
    if (args.length != 2)
      CodeGen.perror "wrong number of arguments\n" + opts.to_s
    end
    options.subdir = args[0]
    options.dirname = File.join "src", options.subdir
    if options.classname
      options.basename = File.join(options.dirname,
                                   File.basename(options.classname) +
                                     "_" + args[1])
    else
      options.basename = File.join options.dirname, args[1]
    end
    return options
  end
end

def create_rep_file(options)
  repfile = File.new(options.repname, 'w')
  repfile << CodeGen.create_license
  guard_macro = File.basename(options.repname).gsub(/\./, '_').upcase
  fkt_pref = options.fkt_pref
  classN = options.classN
  subdir = options.subdir
  content = ERB.new($repfile).result(binding)
  repfile << ERB.new($hwrapper).result(binding)
  repfile.close
end

def fill_out_files(options)
  repfile = File.readlines(options.repname)

  options.functions =
    CodeGen.extract_functions_from_rep(options.classN, repfile)

  repfile.each_with_index do |line, idx|
    next unless line.match /XX\w+XX/
    if line.match /XXfunctionNamesXX/
      repfile[idx] = []
      options.functions.each do |func, parameters|
        repfile[idx].push "  #{options.classN}#{func} #{
          func.gsub(/(^|[a-z])([A-Z])/) do
            replace = ""
            replace += $1 + '_' if $1.length > 0
            replace += $2.downcase
          end};\n"
      end
    end
    if line.match /XXfunctionPtrXX/
      repfile[idx] = []
      indent = repfile[idx-1][/^.*\(/].length
      options.functions.each do |func, parameters|
        repfile[idx].push "#{' ' * indent}#{options.classN}#{func},\n"
      end
      repfile[idx][-1] = repfile[idx][-1].chomp.chop + ");\n"
    end
  end
  repfile.flatten!
  newfile = File.new(options.repname, 'w')
  newfile << repfile.join
  newfile.close
end

def create_interface_header(options)
  options.interfacename = options.basename + ".h"
  if File.exist?(options.interfacename) and not options.force
    CodeGen.perror "file #{options.interfacename} already exists\n" +
           "use --force to overwrite, with --no-create to reuse\n" + 
           "(#$0 -h for details)"
  end

  i_file = File.new(options.interfacename, 'w')
  i_file << CodeGen.create_license

  classN = options.classN
  fkt_pref = options.fkt_pref
  functions = options.functions
  max_type_len = 0
  max_type_len = "#{classN} *".length if options.refcount
  functions.each do |func, parameters|
    max_type_len = parameters[0].length if max_type_len < parameters[0].length
  end

  content = ERB.new($interface_header).result(binding)
  guard_macro = File.basename(options.interfacename).gsub(/\./, '_').upcase
  i_file << ERB.new($hwrapper).result(binding)
end

def create_interface_code(options)
  options.interfacename = options.basename + ".c"
  if File.exist?(options.interfacename) and not options.force
    CodeGen.perror "file #{options.interfacename} already exists\n" +
           "use --force to overwrite, with --no-create to reuse\n" +
           "(#$0 -h for details)"
  end
  i_code_file = File.new(options.interfacename, 'w')
  i_code_file << CodeGen.create_license

  classN = options.classN
  fkt_pref = options.fkt_pref
  functions = options.functions
  interface_funcs = ""
  functions.each do |func, paras|
    next if func == 'DeleteFunc'
    type = paras.shift
    interface_funcs += ERB.new($interface_func).result(binding)
  end
  create_func = ERB.new($create_func).result(binding)
  ref_func = ERB.new($ref_func).result(binding)
  cast_func = ERB.new($cast_func).result(binding)
  delete_func = ERB.new($delete_func).result(binding)
  class_new_fkt = ERB.new($class_new_fkt).result(binding)
  subdir = options.subdir
  i_code_file << ERB.new($interface_file).result(binding)
  i_code_file.close
end

begin
  options = ClassGenOpts.parse(ARGV);
rescue OptionParser::InvalidOption => e
  CodeGen.perror e.message
end

options.repname = "#{options.basename}_rep.h"
options.fkt_pref = File.basename(options.repname, '_rep.h')
options.classN = CodeGen.filename2classname(options.basename)
if File.exist?(options.repname)
  if not options.force
    CodeGen.perror "file #{options.repname} already exists\n" +
           "use --force to overwrite, with --no-create to reuse\n" + 
           "(#$0 -h for details)"
  else
    if options.create
      create_rep_file(options)
    end
  end
else
  create_rep_file(options)
end
system("#{options.editor} #{options.repname}")
if options.create
  fill_out_files(options)
else
  repfile = File.readlines(options.repname)
  options.functions =
    CodeGen.extract_functions_from_rep(options.classN, repfile)
end
create_interface_header(options)
create_interface_code(options)
