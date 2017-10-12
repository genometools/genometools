#!/usr/bin/env ruby

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <indexname> <inputfile>"
  exit 1
end

indexname=ARGV[0]
inputfile=ARGV[1]
puts "env -i bin/gt encseq encode -indexname #{indexname}.unsorted #{inputfile}"
puts "TMPFILE=`mktemp TMP.XXXXXX` || exit 1"
puts "env -i bin/gt seqorder -sortlength #{indexname}.unsorted > \${TMPFILE}"
puts "env -i bin/gt encseq encode -indexname #{indexname} \${TMPFILE}"
puts "rm -f \${TMPFILE} #{indexname}.unsorted.*"
