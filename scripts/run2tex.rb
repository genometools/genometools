#!/usr/bin/env ruby

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <runfile>"
  exit 1
end

def makehashkey(filename,parms)
  return filename + ' ' + parms
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
      if filename == ''
        STDERR.puts "line #{line}: filename is undefined"
        exit 1
      end
      if parms == ''
        STDERR.puts "line #{line}: parms is undefined"
        exit 1
      end
      timeresult=t[1].to_f
      runtimes[makehashkey(filename,parms)] = timeresult
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

def makefastestlist(ftab,ptab,rtimes)
  fastestlist = Hash.new()
  ftab.each_key do |filename|
    timelist = []
    ptab.each_key do |parms|
      t = rtimes[makehashkey(filename,parms)]
      if t
        timelist.push t
      end
    end
    fastestlist[filename] = timelist.min
  end
  return fastestlist
end

numoffiles=filenametab.length
puts "\\begin{tabular}{|l|*{#{numoffiles}}{r|}}\\hline"
puts "parameter & \\multicolumn{#{numoffiles}}{c|}{files}"
puthline
orderedkeys(filenametab).each do |filename|
  print "&#{filename}"
end
puthline
fastestlist=makefastestlist(filenametab,parmstab,runtimes)

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
  end
  puthline
end

print <<'TEXT'
\end{tabular}
\end{center}
\end{document}
TEXT
