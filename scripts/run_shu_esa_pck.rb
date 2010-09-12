#!/usr/bin/ruby
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'evalshurun'
require 'fileutils'

$bin=ENV['GTDIR'] + "/bin/"

def shu_pck(files, param)
  system("#{$bin}gt       " +
           "packedindex mkindex " +
           "-db #{files}        " +
           "-dna                " +
           "-dir rev            " +
           "-ssp                " +
           "-bsize 8            " +
           "-sprank             " +
           "-pl                 " +
           "-sat uchar          " +
           "-indexname pck")
  system("#{$bin}gt genomediff #{param} -pck pck")
end

def shu_esa(files, param)
    system"#{$bin}gt suffixerator -db #{files} -indexname esa " + 
             "-dna -suf -tis -lcp -ssp"
    system"#{$bin}gt genomediff #{param} -esa esa"
end

def reverse_and_concat(file)
  tf = Tempfile.new("gt_rev")
  tmpfile = tf.path
  tf.close
  dirname = File.dirname(file)
  `$GTDIR/bin/gt            \
   convertseq -o #{tmpfile} \
   -force -r #{file}`
  basename = File.basename(file, ".*")
  newfile = File.join(dirname, basename + "_plus_rev.fas")
  FileUtils.cp(file, newfile)
  file = newfile
  File.open(file, 'a') {|fp|
    File.open(tmpfile, 'r') {|tf|
      while line = tf.gets
        fp.puts line
      end
    }
  }
  exit 1 unless File.exist?(newfile)
  return newfile
end

files = ""
ARGV.each do |file|
  unless File.exist?(file)
    puts "non existing file " + file
  else
    files += reverse_and_concat(file) + " "
  end
end

puts "***ESA***"
shu_esa(files, "-v")
system "gprof $GTDIR/bin/gt > esa_full.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK***"
shu_pck(files, "-v")
system "gprof $GTDIR/bin/gt > pck_full.prof"
FileUtils.rm_f("gmon.out")
puts "***ESA***"
shu_esa(files, "-v -shulen")
system "gprof $GTDIR/bin/gt > esa_shulen.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK***"
shu_pck(files, "-v -shulen")
system "gprof $GTDIR/bin/gt > pck_shulen.prof"
FileUtils.rm_f("gmon.out")
puts "***PCK TRAVERSE***"
shu_pck(files, "-v -traverse")
system "gprof $GTDIR/bin/gt > pck_traverse.prof"
FileUtils.rm_f("gmon.out")
