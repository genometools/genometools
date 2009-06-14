#!/usr/bin/env ruby

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <runfile>"
  exit 1
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

fname=ARGV[0]

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
      end
    end
  end
end

def puthline
  puts "\\\\\\hline"
end

print <<'TEXT'
\documentclass[12pt]{article}
\begin{document}
\begin{center}
TEXT

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

numoffiles=filenametab.length
puts "\\begin{tabular}{|l|*{#{2*numoffiles}}{r|}}\\hline"
puts "parameter & \\multicolumn{#{2*numoffiles}}{c|}{files}"
puthline
orderedkeys(filenametab).each do |filename|
  print "&\\multicolumn{2}{|c|}{#{filename}}"
end
puthline
fastestlist=makebestlist(filenametab,parmstab,runtimes)
spaceefflist=makebestlist(filenametab,parmstab,spacereq)

orderedkeys(parmstab).each do |parms|
  print "\\texttt{#{parms}}"
  orderedkeys(filenametab).each do |filename|
    t = runtimes[makehashkey(filename,parms)]
    print "&"
    if t
      if t == fastestlist[filename]
        printf("\\textbf{%.2f}",t)
      else
        printf("%.2f",t)
      end
    end
    print "&"
    s = spacereq[makehashkey(filename,parms)]
    if s
      if s == spaceefflist[filename]
        printf("\\textbf{%.2f}",s)
      else
        printf("%.2f",s)
      end
    end
  end
  puthline
end

print <<'TEXT'
\end{tabular}
\end{center}
\end{document}
TEXT
