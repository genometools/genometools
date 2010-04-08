#!/usr/bin/env ruby

funcs = `nm lib/libgtunstable.so | egrep '^[^ ]{8} T [^g][^t][^_]' | egrep -v '^[^ ]{8} T (XML|md5|zError|zcalloc|zcfree|zlibCompileFlags|zlibVersion|lua|inflate|deflate|crc32|fsetkey|fencrypt|compress|adler32|_fini|_init|Xml)'`.split(/\n/)

funcs = funcs.collect { |f| f.split[2]}.sort{|a,b| b.length <=> a.length}

puts "#{funcs.length} functions to process"

funcs.each do |f|
  (funcs - [f]).each do |f2|
    if f.include?(f2) then
      puts "function list is not prefix-free: #{f} <-> #{f2}"
      #exit
    end
  end
end

funcs.each do |funcname|
  puts funcname
  files = `fgrep -Rl #{funcname} src/*`
  files.each_line do |file|
    fullpath = File.join(`pwd`.chomp,file).chomp
    puts fullpath
    buffer = ''
    File.open(fullpath) do |infile|
      buffer = infile.read
    end
    buffer.gsub!(/(^| |\(|\)|,)(([*!]+)?)#{funcname}([^;])/, "\\1\\2gt_#{funcname}\\4")
    File.open(fullpath, 'w') do |outfile|
      outfile.write(buffer)
    end
  end
end
