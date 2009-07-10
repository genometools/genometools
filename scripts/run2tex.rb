#!/usr/bin/env ruby

def usage()
  STDERR.puts "Usage: #{$0} space|time <runfile>"
  exit 1
end

def puthline
  puts "\\\\\\hline"
end

def makehashkey(filename,parms)
  return filename + ' ' + parms
end

def checkvalues(filename,parms)
  if filename == ''
    STDERR.puts "filename is undefined"
    exit 1
  end
  if parms == ''
    STDERR.puts "parms is undefined"
    exit 1
  end
end

def orderedkeys(ht)
  return ht.keys.sort
end

def makebestlist(ftab,ptab,tab)
  bestlist = Hash.new()
  ftab.each_key do |filename|
    tl = []
    ptab.each_key do |parms|
      t = tab[makehashkey(filename,parms)]
      if t
        tl.push t
      end
    end
    bestlist[filename] = tl.min
  end
  return bestlist
end

def maketab(filenametab,parmstab,datestring,valuehash)
  numoffiles=filenametab.length
  puts "\\begin{sidewaystable}"
  puts "\\begin{center}"
  puts "\\begin{small}"
  puts "\\begin{tabular}{|l|*{#{numoffiles}}{r|}}\\hline"
  puts "parameter & \\multicolumn{#{numoffiles}}{c|}{files}"
  puthline
  orderedkeys(filenametab).each do |filename|
    print "&#{filename}"
  end
  puthline
  bestlist=makebestlist(filenametab,parmstab,valuehash)
  orderedkeys(parmstab).each do |parms|
    print "\\texttt{#{parms}}"
    orderedkeys(filenametab).each do |filename|
	print "&"
	v = valuehash[makehashkey(filename,parms)]
	if v
	  if v == bestlist[filename]
	    printf("\\textbf{%.2f}",v)
	  else
	    printf("%.2f",v)
	  end
	end
    end
    puthline
  end
  print <<'TEXT'
  \end{tabular}
  \end{small}
  \end{center}
TEXT
  puts "\\caption{run at #{datestring}}"
  puts "\\end{sidewaystable}"
end

if ARGV.length != 2
  usage()
end

showtime=false
showspace=false
if ARGV[0] == 'time'
  settime=true
elsif ARGV[0] == 'space'
  setspace=true
elsif ARGV[0] == 'time+space'
  settime=true
  setspace=true
else
  usage()
end
fname=ARGV[1]

begin
  f = File.open(fname,"r")
rescue => err
 STDERR.print "#{$0}: cannot open file \"#{fname}\": #{err}\n"
 exit 1
end

filenametab = Hash.new()
parmstab = Hash.new()
runtimes = Hash.new()
spacereq = Hash.new()
datestring=''

parms=''
filename=''
f.each_line do |line|
  m = line.match(/# RUN ([a-zA-Z0-9]*) (.*)/)
  if m
    filename = m[1]
    filenametab[filename] = true
    parms = m[2]
    if parms == ''
      parms='nothing'
    end
    parmstab[parms] = true
  else
    t = line.match(/# TIME overall ([0-9\.]*)/)
    if t
      checkvalues(filename,parms)
      timeresult=t[1].to_f
      runtimes[makehashkey(filename,parms)] = timeresult
    else
      s = line.match(/# space peak in megabytes: ([0-9\.]*)/)
      if s
        checkvalues(filename,parms)
        space=s[1]
        spacereq[makehashkey(filename,parms)] = space.to_f
      else
        d = line.match(/# DATE ([0-9\-\:]*)/)
        if d
          datestring=d[1]
        end
      end
    end
  end
end

print <<'TEXT'
\documentclass[11pt]{article}
\usepackage{rotating}
\usepackage{a4wide}
\begin{document}
TEXT

if settime
  maketab(filenametab,parmstab,datestring,runtimes)
end

if setspace
  maketab(filenametab,parmstab,datestring,spacereq)
end

print <<'TEXT'
\end{document}
TEXT
