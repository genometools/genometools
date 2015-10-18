#!/usr/bin/env ruby

require 'optparse'
require 'ostruct'
require 'set'

def listdirectory(directory)
  # prepare regexp for entries to ignore: saves time for repeated use
  ignore_dirs = Regexp.compile(/^\.\.?$/)  # match . or ..
  # Iterate over items in directory
  stack = Array.new()
  stack.push(directory)
  filelist = Array.new()
  while not stack.empty? do
    currentdir = stack.pop
    Dir.foreach(currentdir) do |entry|
      if not entry.match(ignore_dirs)
      # directory entry is a regular file
        if File.stat("#{currentdir}/#{entry}").file?
          filelist.push("#{currentdir}/#{entry}")
        # directory entry is a subdirectory
        elsif File.stat("#{currentdir}/#{entry}").directory?
          stack.push("#{currentdir}/#{entry}")
        end
      end
    end
  end
  return filelist
end

def showresult(prog,minidentity,found,fails)
  print "minid=#{minidentity} #{prog} #{found} of #{found+fails}, "
  printf("miss #{fails}, sensitivity %.2f%%\n",
         100.0 * found.to_f/(found+fails).to_f)
end

def matchinfile(filename)
  if open(filename).grep(/^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/).length > 0
    return true
  else
    return false
  end
end

def indentify_num(regexp,s)
  if m = s.match(regexp)
    return m[1].to_i
  end
  STDERR.puts "#{$0}: cannot identify #{regexp} in #{s}"
  exit 1
end

def makesystemcall(argstring)
  if not system(argstring)
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1
  end
end

def makedirname(dir,minidentity)
  return "#{dir}/minid-#{minidentity}"
end

def makefilename(dir,minidentity,pref,length,seqnum)
  return "#{makedirname(dir,minidentity)}/#{pref}-#{length}-#{seqnum}"
end

def runseedextend(inputdir,targetdir,seedlength,minidentity,length,seqnum)
  gtcall = "env -i bin/gt"
  inputfile = makefilename(inputdir,minidentity,"rand",length,seqnum)
  destfile = makefilename(targetdir,minidentity,"gtout",length,seqnum)
  destdir = File.dirname(destfile)
  makesystemcall("mkdir -p #{destdir}")
  makesystemcall("#{gtcall} encseq encode -des no -sds no -md5 no " +
                 "-indexname #{inputfile} #{inputfile}.fas")
  makesystemcall("#{gtcall} seed_extend -t 21 -l #{length} " +
		 "-seedlength #{seedlength} -minidentity #{minidentity} " +
		 "-seed-display -bias-parameters -extendgreedy " +
		 "-overlappingseeds -ii #{inputfile} > " +
		 "#{destfile}.txt")
  succ = matchinfile("#{destfile}.txt")
  ["esq","ssp"].each do |suffix|
    filename = "#{inputfile}.#{suffix}"
    if File.exists?(filename)
      File.delete(filename)
    end
  end
  return succ
end

def rundaligner(inputdir,targetdir,seedlength,minidentity,length,seqnum)
  inputfile = makefilename(inputdir,minidentity,"rand",length,seqnum)
  destfile = makefilename(targetdir,minidentity,"daout",length,seqnum)
  destdir = File.dirname(destfile)
  makesystemcall("mkdir -p #{destdir}")
  makesystemcall("scripts/convert2myersformat.rb #{inputfile}.fas " +
		 "> #{destfile}.fasta")
  if ENV.has_key?("PACKAGES")
    myersprog = ENV["PACKAGES"]
  else
    myersprog = ".."
  end
  makesystemcall("#{myersprog}/DAZZ_DB/fasta2DB #{destfile}.db #{destfile}.fasta")
  makesystemcall("#{myersprog}/DALIGNER/daligner -t21 -I -A -Y " +
                 "-e0.#{minidentity} -l#{length} -k#{seedlength} " +
                 "#{destfile}.db #{destfile}.db > #{destfile}.txt")
  File.delete("#{destfile}.db")
  return matchinfile("#{destfile}.txt")
end

def parseargs(argv)
  options = OpenStruct.new
  options.runse = nil
  options.runda = nil
  options.num_tests = nil
  options.minid_max = nil
  options.minid_min = nil
  options.inputdir = ENV["HOME"] + "/dalign-files"
  options.targetdir = "dalign"
  options.first = 0
  opts = OptionParser.new
  indent = " " * 37
  opts.on("-g","--generate-seq STRING",
          "generate-sequences argument:" +
          "\n#{indent}argument: num_test,minid_min,minid_max") do |x|
    runseqgen = x.split(/,/).map {|x| x.to_i}
    if runseqgen.length != 3
      STDERR.puts "#{$0}: need three integer arguments for option -q"
      exit 1
    end
    options.num_tests = runseqgen[0]
    options.minid_min = runseqgen[1]
    options.minid_max = runseqgen[2]
  end
  opts.on("-s","--seed_extend","run gt seed_extend") do |x|
    options.runse = true
  end
  opts.on("-d","--daligner","run daligner") do |x|
    options.runda = true
  end
  opts.on("-i","--inputdir STRING","specify input directory" +
          "\n#{indent}(default: #{options.inputdir})") do |x|
    options.inputdir = x
  end
  opts.on("-t","--targetdir STRING","specify target directory" +
               "\n#{indent}(default: #{options.targetdir}") do |x|
    options.targetdir = x
  end
  opts.on("-f","--first NUMBER",
          "specify number of sequence used for evaluation") do |x|
    options.first = x.to_i
    if options.first <= 0
      STDERR.puts "#{$0}: argument of option -f must be positive"
      exit 1
    end
  end
  rest = opts.parse(argv)
  if rest.length != 0 or (options.num_tests.nil? and not options.runda and 
                          not options.runse)
    STDERR.puts "#{opts}"
    exit 1
  end
  return options
end

options = parseargs(ARGV)
seedlength = 14

if not options.minid_min.nil?
  if options.num_tests < 20 then
    lengthtab = [*90..110].sample(options.num_tests)
  else
    lengthtab = options.num_tests.times.map{ 90 + Random.rand(21) }
  end
  options.minid_min.upto(options.minid_max).each do |minidentity|
    filedir = makedirname(options.targetdir,minidentity)
    makesystemcall("mkdir -p #{filedir}")
    lengthtab.each_with_index do |length,seqnum|
      # generate random sequences
      outfile = "#{filedir}/rand-#{length}-#{seqnum}"
      makesystemcall("ruby scripts/gen-randseq.rb " +
                     "--minidentity #{minidentity} " +
		     "--seedlength #{seedlength} --length #{length+10} -c 35 " +
		     "--mode seeded 1> #{outfile}.fas 2>/dev/null")
    end
  end
end
if options.runse or options.runda
  gtfound = Array.new(100) {0}
  gtfail = Array.new(100) {0}
  dafound = Array.new(100) {0}
  dafail = Array.new(100) {0}
  minidset = Set.new()
  listdirectory(options.inputdir).each do |filename|
    if filename.match(/\.fas/)
      seqnum = indentify_num(filename,Regexp.compile(/rand-\d+-(\d+)\./))
      if options.first == 0 or seqnum < options.first
        minidentity = indentify_num(filename,Regexp.compile(/minid-(\d+)\//))
        minidset.add(minidentity)
        length = indentify_num(filename,Regexp.compile(/rand-(\d+)-/))
        if options.runse
          if runseedextend(options.inputdir,options.targetdir,seedlength,
                           minidentity,length,seqnum)
            gtfound[minidentity] += 1
          else
            gtfail[minidentity] += 1
          end
        end
        if options.runda
          if rundaligner(options.inputdir,options.targetdir,seedlength,
                         minidentity,length,seqnum)
            dafound[minidentity] += 1
          else
            dafail[minidentity] += 1
          end
        end
      end
    end
  end
  minidset.sort.each do |minidentity|
    if options.runse
      showresult("seed_extend",minidentity,gtfound[minidentity],
                                           gtfail[minidentity])
    end
    makesystemcall("rm -f daout-*.daout*.las")
    if options.runda
      showresult("daligner",minidentity,dafound[minidentity],
                                        dafail[minidentity])
    end
  end
end
