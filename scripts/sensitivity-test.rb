#!/usr/bin/env ruby

require 'optparse'
require 'ostruct'
require 'set'

def indentify_num(s,regexp)
  if m = s.match(regexp)
    return m[1].to_i
  end
  STDERR.puts "#{$0}: cannot identify #{regexp} in #{s}"
  exit 1
end

def fromfilename2keys(filename)
  minidentity = indentify_num(filename,Regexp.compile(/minid-(\d+)\//))
  length = indentify_num(filename,Regexp.compile(/-(\d+)-/))
  seqnum = indentify_num(filename,Regexp.compile(/-\d+-(\d+)\./))
  return minidentity, length, seqnum
end

def makesystemcall(argstring,withecho = false)
  if not system(argstring)
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1
  elsif withecho
    puts "# #{argstring}"
  end
end

def matchinfile(filename)
  if open(filename).grep(/^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/).length > 0
    return true
  else
    return false
  end
end

def makedirname(dir,minidentity)
  return "#{dir}/minid-#{minidentity}"
end

def makefilename(dir,minidentity,pref,length,seqnum)
  return "#{makedirname(dir,minidentity)}/#{pref}-#{length}-#{seqnum}"
end

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

def gtcall ()
  return "env -i bin/gt"
end

def callseedextend(indexname,inputfile,destfile,minidentity,length,seqnum,
                   seedlength,weakends,withalignment,withecho = false)
  makesystemcall("#{gtcall()} encseq encode -des no -sds no -md5 no " +
                 "-indexname #{indexname} #{inputfile}.fas",withecho)
  makesystemcall("#{gtcall()} seed_extend -t 21 -l #{length} " +
		 "-seedlength #{seedlength} -minidentity #{minidentity} " +
		 "-seed-display -bias-parameters -extendgreedy " +
		 "-overlappingseeds -ii #{indexname}" +
                 (if withalignment then " -a" else "" end) +
                 (if weakends then " -weakends" else "" end) +
		 (if destfile.empty? then "" else " > #{destfile}.txt" end),
                 withecho)
  succ = true
  if destfile != ""
    succ = matchinfile("#{destfile}.txt")
    ["esq","ssp"].each do |suffix|
      filename = "#{indexname}.#{suffix}"
      if File.exists?(filename)
        File.delete(filename)
      end
    end
  end
  return succ
end

def runseedextend(inputdir,targetdir,weakends,withalignment,seedlength,
                  minidentity,length,seqnum)
  inputfile = makefilename(inputdir,minidentity,"rand",length,seqnum)
  destfile = makefilename(targetdir,minidentity,"gtout",length,seqnum)
  destdir = File.dirname(destfile)
  makesystemcall("mkdir -p #{destdir}")
  return callseedextend(inputfile,inputfile,destfile,minidentity,length,seqnum,
                        seedlength,weakends,withalignment)
end

def rundaligner(inputdir,targetdir,seedlength,minidentity,length,seqnum,
                tofile = true)
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
  withecho = if tofile then false else true end
  if tofile
    outputfile = "> #{destfile}.txt"
  else
    outputfile = ""
    makesystemcall("#{gtcall()} encseq encode -des no -sds no -md5 no " +
                   "-indexname #{destfile} #{destfile}.fasta",withecho)
  end
  if withecho
    puts "# daligner result:"
  end
  makesystemcall("#{myersprog}/DALIGNER/daligner -t21 -I -A -Y " +
                 "-e0.#{minidentity} -l#{length} -k#{seedlength} " +
                 "#{destfile}.db #{destfile}.db #{outputfile}",withecho)
  if not tofile
    makesystemcall("#{gtcall()} dev show_seedext -f #{destfile}.txt -a " +
                   " -polished-ends",
                   withecho)
  end
  File.delete("#{destfile}.db")
  return matchinfile("#{destfile}.txt")
end

def rerun_seedextend(options,seedlength)
  filename = options.rerun
  minidentity, length, seqnum = fromfilename2keys(filename)
  inputfiledir = options.inputdir
  puts # "minid=#{minidentity}, length=#{minidentity}, seqnum=#{seqnum}"
  inputfile = makefilename(inputfiledir,minidentity,"rand",length,seqnum)
  puts "# inputfile=#{inputfile}.fas"
  indexname = "sfx-#{length}-#{seqnum}"
  callseedextend(indexname,inputfile,"",minidentity,length,seqnum,
                 seedlength,options.weakends,true,true)
  rundaligner(options.inputdir,options.targetdir,seedlength,
              minidentity,length,seqnum,false)
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
  options.rerun = nil
  options.first = 0
  options.weakends = false
  opts = OptionParser.new
  indent = " " * 37
  opts.banner = "Usage: #{$0} [options] \nFirstly, generate sequences using "+
  "the -g option specifying the number of samples\nand the range of "+
  "minidentity of the sequences. With -t you should set the target\ndirectory "+
  "to #{options.inputdir}.\nAfter that, execute this script again to run "+
  "gt_seed_extend (-s) and/or daligner\n(-d) on the testdata. For a testdata "+
  "subset specify a number with opton -f.\nThe DALIGNER and DAZZ_DB "+
  "directories required by -d are expected in ../ but you\ncan specify a "+
  "different location with environment variable ${PACKAGES}."

  opts.on("-g","--generate-seq STRING",
          "generate-sequences argument:" +
          "\n#{indent}argument: num_test,minid_min,minid_max") do |x|
    runseqgen = x.split(/,/).map {|x| x.to_i}
    if runseqgen.length != 3
      STDERR.puts "#{$0}: need three integer arguments for option -g"
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
  opts.on("-w","--weakends","use option -weakends of seed_extend") do |x|
    options.weakends = true
  end
  opts.on("-i","--inputdir STRING","specify input directory" +
          "\n#{indent}(default: #{options.inputdir})") do |x|
    options.inputdir = x
  end
  opts.on("-t","--targetdir STRING","specify target directory" +
               "\n#{indent}(default: #{options.targetdir})") do |x|
    options.targetdir = x
  end
  opts.on("-r","--rerun STRING","rerun seed-extend on specific file") do |x|
    options.rerun = x
  end
  opts.on("-f","--first NUMBER",
          "specify number of sequences used for\n#{indent}evaluation") do |x|
    options.first = x.to_i
    if options.first <= 0
      STDERR.puts "#{$0}: argument of option -f must be positive"
      exit 1
    end
  end
  opts.on("-h", "--help", "print this help message") do
    puts opts
    exit
  end
  rest = opts.parse(argv)
  if rest.length != 0 or (options.num_tests.nil? and not options.runda and
                          not options.runse and options.rerun.nil?)
    STDERR.puts "#{opts}"
    exit 1
  end
  if options.weakends and not options.seed_extend
    STDERR.puts "option -w/--weakends requires option -s/--seed_extend"
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
  makesystemcall("mkdir -p #{options.targetdir}")
  listdirectory(options.inputdir).each do |filename|
    if filename.match(/\.fas/)
      minidentity, length, seqnum = fromfilename2keys(filename)
      if options.first == 0 or seqnum < options.first
        minidset.add(minidentity)
        if options.runse
          if runseedextend(options.inputdir,options.targetdir,options.weakends,
                           false,seedlength,minidentity,length,seqnum)
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
  puts "# " + ARGV.join(" ")
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
if not options.rerun.nil?
  rerun_seedextend(options,seedlength)
end
