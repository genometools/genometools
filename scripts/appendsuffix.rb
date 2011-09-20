#!/usr/bin/ruby

#
# Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
# Copyright (c) 2011 Center for Bioinformatics, University of Hamburg
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

usage = <<-end_usage
Usage:
  #$0 <templatefilename> <suffix> [<identifier>+]

Find macro names, function names and typedefs in a C file (templatefilename)
and append "_" followed by the specified suffix to these identifiers and to
further identifiers specified as optional arguments.

end_usage

if ARGV.size < 2 || !File.exists?(ARGV[0])
  puts usage
  exit
end

class IdentifiersCollection

  attr_reader :macros, :functions, :types

  def initialize(c_code)
    uncommented = delcomments(c_code)
    toplevel = delblocks(uncommented)
    @macros = findmacronames(toplevel)
    notpreproc_toplevel = delpreproc(toplevel)
    @functions = findfunctionnames(notpreproc_toplevel)
    @types = findtypes(notpreproc_toplevel)
  end

  private

  def delcomments(input)
    output = ""
    prevc = nil
    state = :code
    input.each_byte do |c|
      if state == :comment
        state = :code if c == ?/ and prevc == ?*
      else
        if c == ?* and prevc == ?/
          state = :comment
          output.chop!
        else
          output << c
        end
      end
      prevc = c
    end
    return output
  end

  def delblocks(uncommented)
    output = ""
    level = 0
    uncommented.each_byte do |c|
      level += 1 if c == ?{
      output << c if level == 0
      level -= 1 if c == ?}
    end
    return output
  end

  def delpreproc(input)
    state = :notpreproc
    prevc = nil
    output = ""
    input.each_byte do |c|
      if state == :preproc
        state = :notpreproc if c == ?\n and prevc != ?\\
      else
        if c == ?#
          state = :preproc
        else
          output << c
        end
      end
      prevc = c
    end
    return output
  end

  def findmacronames(toplevel)
    toplevel.scan(/#\s*define\s(\w+)/).flatten
  end

  def findfunctionnames(notpreproc_toplevel)
    notpreproc_toplevel.scan(/(\w+)\s*\([^\*].*?\)/m).flatten
  end

  def findtypes(notpreproc_toplevel)
    notpreproc_toplevel.scan(/typedef\s+[^;]*\s+(\w+)\s*;/m).flatten
  end

end

class CCode

  def initialize(c_code)
    @code = c_code
    @ids = IdentifiersCollection.new(c_code)
  end

  def append_suffix!(suffix, more_identifiers = [])
    ids = more_identifiers
    [:macros, :functions, :types].each do |idtype|
      ids += @ids.send(idtype)
    end
    ids.each do |id|
      @code.gsub!(/(\W+)(#{id})(\W+)/m) {"#$1#$2_#{suffix}#$3"}
    end
  end

  def to_s
    @code
  end

end

input = IO.read(ARGV.shift)
suffix = ARGV.shift
more_identifiers = ARGV

code = CCode.new(input)
code.append_suffix!(suffix, more_identifiers)
print code
