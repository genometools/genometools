#!/usr/bin/ruby
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'evalshurun'

$bin=ENV['GTDIR'] + "/bin/"

def shu_pck(files, param)
  system("#{$bin}gt       " +
           "packedindex mkindex " +
           "-db #{files}        " +
           "-dna                " +
           "-dir rev            " +
           "-ssp                " +
           "-dc 64              " +
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
    files += file + " "
  end
end

puts "***ESA***"
shu_esa(files, "-v")
puts "***PCK***"
shu_pck(files, "-v")
puts "***ESA***"
shu_esa(files, "-v -shulen")
puts "***PCK***"
shu_pck(files, "-v -shulen")
puts "***PCK TRAVERSE***"
shu_pck(files, "-v -traverse")
