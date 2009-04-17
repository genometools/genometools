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

def each_seq(filename)
  File.open(filename) do |file|
    (farr = file.read.split(">")).each do |fitem|
      seq = fitem.split("\n")
      fline = seq.shift
      next if fline.nil?
      sequence = seq.collect{|l| l.chomp}.join
      yield fline, sequence
    end
  end
end

if ARGV.length != 1 then
  STDERR.puts "Usage: #{$0} <FASTA file>"
  exit 1
end

j=0
begin
  each_seq(ARGV[0]) do |header, seq|
    # break up headers into words and create description lines not longer than
    # 70 characters of length
    hlines = [[]]
    llen = 0
    lineno = 0
    header.chomp.split(/\s+/).each do |word|
      if llen + word.length + 1 > 70
        hlines.push([])
        lineno += 1
        llen = 0
      end
      hlines[lineno].push(word)
      llen += word.length + 1
    end
    i=1
    puts "ID   sequence#{j}"
    puts "XX"
    hlines.each do |hline|
      puts "DE   #{hline.join(" ")}"
    end
    puts "XX"
    puts "SQ"
    print "     "
    seq.each_byte do |c|
      putc c
      if i % 10 == 0 then
        print " "
      end
      if i % 60 == 0 then
        printf("%9s\n     ", i)
      end
      i += 1
    end
    printf(' '*(80-i%60-(i%60)/10-13) + "%9d\n", i-1)
    puts "//\n\n"
    j+=1
  end
rescue StandardError => msg
  STDERR.puts "Error: #{msg}"
  exit 1
end
