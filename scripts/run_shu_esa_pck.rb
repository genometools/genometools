#!/usr/bin/ruby
#
# Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'genomediff'
require 'fileutils'

files = ""
ARGV.each do |file|
  unless File.exist?(file)
    puts "non existing file " + file
  else
    files += Genomediff.reverse_and_concat(file) + " "
  end
end

puts "***ESA***"
Genomediff.esa_index(files, 1, "esa")
Genomediff.esa_genomediff("esa"," -v")
system "gprof $GTDIR/bin/gt > esa_full.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK***"
Genomediff.pck_index(files, 8, 1, "pck")
Genomediff.pck_genomediff("pck"," -v")
system "gprof $GTDIR/bin/gt > pck_full.prof"
FileUtils.rm_f("gmon.out")
puts "***ESA***"
Genomediff.esa_genomediff("esa"," -v -shulen")
system "gprof $GTDIR/bin/gt > esa_shulen.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK***"
Genomediff.pck_genomediff("pck"," -v -shulen")
system "gprof $GTDIR/bin/gt > pck_shulen.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK TRAVERSE***"
Genomediff.pck_genomediff("pck"," -v -traverse")
system "gprof $GTDIR/bin/gt > pck_traverse.prof"
FileUtils.rm_f("gmon.out")
