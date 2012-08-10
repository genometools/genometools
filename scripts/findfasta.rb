#!/usr/bin/ruby

def listdirectory(directory)
  # prepare regexp for entries to ignore
  # saves time for repeated regexp use, since it stays the same
  ignore_dirs = Regexp.compile(/^\.\.?$/)
  stack = Array.new
  stack.push(directory)
  while not stack.empty?
    d = stack.pop
    Dir.foreach(d) do |entry|
      if not ignore_dirs.match(entry)
        if File.stat("#{d}/#{entry}").file?
          yield "#{d}/#{entry}"
        else
          stack.push("#{d}/#{entry}")
        end
      end
    end
  end
end

def listselected(dirname,excludelist)
  suffixes = ["fastq","fasta","fna","fa","fsa.gz","fsa","FASTA.gz","FASTA"]
  listdirectory(dirname) do |filename|
    suffixes.each do |suffix|
      if filename.match(/\.#{suffix}$/) and
         not excludelist.member?(File.basename(filename))
        yield filename
      end
    end
  end
end

testdata_exclude = ["solid_color_reads.fastq",
                    "test2_wrong_begin.fastq",
                    "test9_uneven_length.fastq",
                    "test7_empty_seq.fastq",
                    "test6_premature_end.fastq",
                    "test4_different_seqlengths.fastq",
                    "test3_different_seqnames.fastq",
                    "corruptpatternfile.fna",
                    "sw100K1.fsa",
                    "sw100K2.fsa"]

listselected("testdata",testdata_exclude) do |filename|
  puts filename
end

if ENV.has_key?("GTTESTDATA")
  gttestdata_exclude = ["trembl-section.fsa.gz"]
  listselected(ENV["GTTESTDATA"],gttestdata_exclude) do |filename|
    puts filename
  end
end
