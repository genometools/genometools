#!/usr/bin/env ruby

GC_code = Struct.new("GC_code",:name,:idnum,:aa,:start)

def gc_code_pretty(gc)
  l = ["  {\"#{gc.name}\"",
       "(unsigned int) #{gc.idnum}",
       "\"#{gc.aa}\"",
       "\"#{gc.start}\"}"]
  return l.join(",\n   ")
end

def transnum_idx_pretty(idx_map,idx)
  if idx_map.has_key?(idx)
    return "  #{idx_map[idx]}U"
  else
    return "  GT_UNDEFTRANSNUM"
  end
end

start_parse = false
current = nil
codelist = Array.new()
idx_map = Hash.new()
maxnum = nil
idx = 0
STDIN.each_line do |line|
  if start_parse
    if m = line.match(/name \"([^\"]*)\"/)
      if not m[1].match(/SGC[0-9]/)
        if not current.nil?
          codelist.push(current)
        end
        current = GC_code.new(m[1],nil,nil,nil)
      end
    elsif m = line.match(/id (\d+)/)
      current.idnum = m[1].to_i
      idx_map[current.idnum] = idx
      idx += 1
      maxnum = current.idnum
    elsif m = line.match(/sncbieaa \"([^\"]*)\"/)
      current.start = m[1]
    elsif m = line.match(/ncbieaa  \"([^\"]*)\"/)
      current.aa = m[1]
    end
  elsif line.match(/^Genetic-code-table/)
    start_parse = true
  end
end
codelist.push(current)

puts "static GtTranslationScheme schemetable[] = {"
puts codelist.map {|gc| gc_code_pretty(gc)}.join(",\n")
puts "};"

puts "\nstatic unsigned int transnum2index[] =\n{"
puts (0..maxnum).to_a.map {|idx| transnum_idx_pretty(idx_map,idx)}.join(",\n")
puts "};"
