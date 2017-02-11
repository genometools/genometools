#!/usr/bin/env ruby

def diag(val)
  return val.bend - val.aend
end

def compare_seeds(val1,val2)
  diff1 = diag(val1)
  diff2 = diag(val2)
  if diff1 < diff2
    return -1
  end
  if diff1 > diff2
    return 1
  end
  if val1[1] < val2[1]
    return -1
  end
  if val1[1] > val2[1]
    return 1
  end
  return 0
end

Maxmatch = Struct.new("Maxmatch",:aend,:bend,:len)

def compare_maxmatches(val1,val2)
  if val1.bend < val2.bend
    return -1
  end
  if val1.bend > val2.bend
    return 1
  end
  if val1.aend < val2.aend
    return -1
  end
  if val1.aend > val2.aend
    return 1
  end
  return 0
end

seeds = Array.new()
seedlength = 15
STDIN.each_line do |line|
  a = line.split(/\s/)
  seeds.push(Maxmatch.new(a[0].to_i,a[1].to_i,seedlength))
end
maxmatches = Array.new()
previous = nil
seeds.sort {|a,b| compare_seeds(a,b)}.each do |val|
  if not previous.nil? 
    if diag(previous) == diag(val) and
       previous.aend + 1 == val.aend and 
       previous.bend + 1 == val.bend
       previous.aend += 1
       previous.bend += 1
       previous.len += 1
    else
      extension = previous.len - seedlength
      if extension < 0
        STDERR.puts "not expected"
        exit 1
      end
      maxmatches.push(Maxmatch.new(previous.aend - extension,
                                   previous.bend - extension,
                                   extension))
      previous = val
    end
  else
    previous = val
  end
end
if not previous.nil?
  extension = previous.len - seedlength
  if extension < 0
    STDERR.puts "not expected"
    exit 1
  end
  maxmatches.push(Maxmatch.new(previous.aend - extension,
                               previous.bend - extension,
                               extension))
end

sumextensions = 0
maxmatches.sort {|a,b| compare_maxmatches(a,b)}.each do |val|
  sumextensions += val.len
  puts "#{val.aend} #{val.bend} #{val.len}"
end
puts "# #{seeds.length} seeds"
puts "# #{maxmatches.length} maxmatches"
if sumextensions + maxmatches.length != seeds.length
  STDERR.puts "not expected"
  exit 1
end
