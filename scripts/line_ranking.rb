#!/usr/bin/env ruby
#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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

require "find"

FILETYPES = %w(c h gen rb py)
IGNOREDIRS = %w(external .git www doc)
IGNOREPATTERNS = [/^\s*\{\s*$/, /^\s*\}\s*$/, /^\s*#/]
COPYRIGHTPATTERN = /Copyright \([cC]\)\s+\d{4}([-, ]*\d{4})?\s+([^<]+)\s+(<)/
C_COMMENT_PATTERN = Regexp.compile('(/\*.*?\*/)', Regexp::MULTILINE)

class Float
  def round_to(x)
    (self * 10**x).round.to_f / 10**x
  end
end

def each_code_file(root)
  Find.find(root) do |path|
    if IGNOREDIRS.include?(File.basename(path.downcase)) then
      STDERR.puts "skipping subtree #{path}"
      Find.prune
    else
      filetype = File.basename(path.downcase).split(".").last
      if FILETYPES.include?(filetype) then
        File.open(path) do |file|
          begin
            yield file
          rescue => msg
            STDERR.printf("#{msg}\n", path)
          end
        end
      end
    end
  end
end

def each_name(str)
  str.scan(COPYRIGHTPATTERN) do |_, name|
    if !name.include?("Center for ") and !name.include?("University of Cali")
      yield name
    end
  end
end

if ARGV.length != 1 then
  STDERR.puts "Usage: #{$0} <source root directory>"
  exit(1)
end

rootdir = ARGV[0]
totallines = 0
users = Hash.new{|hash, key| hash[key] = {:files => 0, :lines => 0}}

each_code_file(rootdir) do |file|
  usersforfile = []
  filecontents = file.readlines
  # find names in source code
  each_name(filecontents.join) do |name|
    usersforfile.push(name)
  end
  # filter out lines with curly braces, one-line comments, ...  
  filecontents.reject! do |line|
    rejected = false
    IGNOREPATTERNS.each do |pat|
     if !rejected and pat.match(line.strip) then
       rejected = true
       break
      end
    end
    rejected
  end
  # do not count C-style multiline comments either
  numlines = filecontents.join.gsub(C_COMMENT_PATTERN, '').split("\n").length
  # credit each author only once per file
  usersforfile.uniq!
  if usersforfile.length == 0 then
    raise "warning: cannot find authors in file %s"
  else
    usersforfile.each do |name|
      users[name][:files] += 1
      # evenly distribute line share between joint authors
      users[name][:lines] += numlines/(usersforfile.length)
    end
    totallines += numlines
  end
end

puts "Code line ranking:"
puts "------------------"
users.sort{|a,b| b[1][:lines] <=> a[1][:lines]}.each do |username, counts|
  puts "#{username}: #{counts[:files]} files, #{counts[:lines]} lines " + \
       "(#{((counts[:lines].to_f/totallines.to_f)*100).round_to(2)}%)"
end
