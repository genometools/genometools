#!/usr/bin/env ruby

def each_tool(binary)
  vals = `#{binary} -help`.split(/\n/)
  loop do
    w = vals.shift
    break if w == "Tools:"
  end
  vals.shift

  tools = []
  loop do
    w = vals.shift
    break if w.nil? or w.strip == ""
    tools << w.split(/\s+/)[0]
  end

  tools.each do |t|
    yield t
  end
end

def get_manpage_for_help(helptext, t)
  fname = "gt-#{t}.1.txt"
  File.open(fname, "w") do |f|
    outlines = helptext.gsub(/-{3,}/,'').split(/\n/)
    short = outlines[1].tr("<>", "[]")
    synopsis = "gt " + outlines[0][7+ARGV[0].length+1..outlines[0].length-1]
    outlines.shift
    outlines.shift
    outlines.shift
    opts = []
    loop do
      w = outlines.shift
      break if w.nil? or w.length  == 0 or !(/[ -]/.match(w[0].chr))
      opts << w
    end

    desc = []
    author = nil
    loop do
      w = outlines.shift
      if (/Report bugs to/).match(w)
        author = w
      end
      break if w.nil? or (/Report bugs to/).match(w)
      desc << w
    end

    if author.nil? then
      loop do
        w = outlines.shift
        if (/Report bugs to/).match(w)
          author = w
        end
        break if w.nil? or (/Report bugs to/).match(w)
      end
    end

    opts = opts.join("\n").gsub("\n-","\n\n-")
    opts.gsub!(/^(-[^ ]+)/, "*\\1*::\n")
    opts.gsub!(/(default: .+$)/, "(\\1)")

    desc.collect! { |l| l.scan(/.{1,79}/)}.join("\n")
    desc = desc.join("\n")
    opts.gsub!(/.*\n-+\n/, "")

    f.puts "GT#{"-#{t.upcase}" if !t.nil?}(1)"
    f.puts "#{if !t.nil? then '=' * t.length end}======"
    f.puts
    f.puts "NAME"
    f.puts "----"
    f.puts "gt#{"-#{t}" if !t.nil?} - #{short}"
    f.puts
    f.puts
    f.puts "SYNOPSIS"
    f.puts "--------"
    f.puts synopsis
    f.puts
    f.puts
    f.puts "OPTIONS"
    f.puts "-------"
    f.puts opts
    f.puts
    f.puts
    if desc.length > 0
      f.puts "DESCRIPTION"
      f.puts "-----------"
      f.puts desc
      f.puts
      f.puts
    end
    if !author.nil?
      f.puts "AUTHOR"
      f.puts "------"
      f.puts author
    end
  end
  `a2x --doctype manpage --format manpage #{fname}`
  if $? == 0 then
    File.unlink(fname)
  end
end

gttext = `#{ARGV[0]} -help`
get_manpage_for_help(gttext, nil)

threads = []

each_tool(ARGV[0]) do |tool|
  threads << Thread.new(tool) do |t|
    text = `#{ARGV[0]} #{t} -help`
    if text.match("Tools:") then
      each_tool("#{ARGV[0]} #{t}") do |t2|
        text2 = `#{ARGV[0]} #{t} #{t2} -help`
        get_manpage_for_help(text2, "#{t}-#{t2}")
      end
    end
    get_manpage_for_help(text, t)
  end
end
threads.each { |aThread|  aThread.join }
