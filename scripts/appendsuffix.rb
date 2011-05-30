#!/usr/bin/env ruby

usage = <<-end_usage
Usage:
  #$0 <templatefilename> <suffix>

Find macros, functions and typedefs in a C file (templatefilename)
and append the given suffix according to the following rules:

  <MACRONAME>_SUFFIX
  <typename>Suffix
  <function_name>_suffix

end_usage

if ARGV.size != 2 || !File.exists?(ARGV[0])
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

  def append_suffix!(suffix)
    suf = {}
    suf[:macros] = "_#{suffix.upcase}"
    suf[:functions] = "_#{suffix.downcase}"
    suf[:types] = suffix.gsub(/^(.)/){$1.upcase}
    [:macros, :functions, :types].each do |idtype|
      @ids.send(idtype).each do |id|
        @code.gsub!(/(\W*)(#{id})(\W*)/) {"#$1#$2#{suf[idtype]}#$3"}
      end
    end
  end

  def to_s
    @code
  end

end

input = IO.read(ARGV.shift)
suffix = ARGV.shift

code = CCode.new(input)
code.append_suffix!(suffix)
print code
