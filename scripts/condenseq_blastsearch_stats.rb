#!/usr/bin/env ruby
# this script compares the blastp and compsearch results, calculating
# sensitivity, precision, tp, fp and fn rates.

class Hit
  include Comparable
  attr_reader :qseqid
  attr_reader :sseqid
  attr_reader :qstart
  attr_reader :qend
  attr_reader :sstart
  attr_reader :send
  attr_reader :length

  def initialize(hitline = "idq ids 0 0 0 0 0 0 0 0")
    parts = hitline.split(" ")
    @qseqid = parts[0]
    @sseqid = parts[1]
    @pident = Float(parts[2])
    @length = Integer(parts[3])
    # if match is on reverse strand, start < end
    # to make overlap calculation more easy, sort
    # strand information is not important here.
    @qstart, @qend = [Integer(parts[4]), Integer(parts[5])].sort
    @sstart, @send = [Integer(parts[6]), Integer(parts[7])].sort
    @evalue = Float(parts[8])
    @bitscore = Float(parts[9])
  end

  def <=> (other)
    diff = @sseqid <=> other.sseqid
    return diff if diff != 0

    diff = @qseqid <=> other.qseqid
    return diff if diff != 0

    diff = @sstart - other.sstart
    return diff if diff != 0

    diff = @send - other.send
    return diff if diff != 0

    diff = @qstart - other.qstart
    return diff if diff != 0

    return @qend - other.qend
  end

  def to_s
    "#@qseqid\t#@sseqid\t#@pident\t#@length\t#@qstart\t#@qend\t#@sstart\t#@send" \
      "\t#@evalue\t#@bitscore"
  end

  def disjuct?(other)
    if (@qseqid != other.qseqid or @sseqid != other.sseqid) or
       (@sstart > other.send or @send < other.sstart) or
       (@qstart > other.qend or @qend < other.qstart)
      return true
    end
    return false
  end

  def overlapping?(other)
    return (not self.disjuct?(other))
  end

end

class HitCollection
  attr_reader :sorted
  include Enumerable
  def initialize(file)
    @sorted = false
    @arr = []
    blast_out = File.new(file, "r")
    blast_out.each_line do |line|
      next if line.match /^#/
      @arr.push Hit.new(line)
    end
    blast_out.close
  end

  def each
    @arr.each do |elem|
      yield elem
    end
  end

  def sort!
    return self if @sorted
    @arr.sort!
    @sorted = true
    return self
  end

  def length
    @arr.length
  end

  def [](idx)
    @arr[idx]
  end

  def &(other)
    new = []
    self.sort!
    other.sort!
    selfidx = 0
    otheridx = 0
    selflen = self.length
    otherlen = other.length

    while selfidx < selflen and otheridx < otherlen
      diff = 0
      if not @arr[selfidx].overlapping?(other[otheridx])
        diff = @arr[selfidx] <=> other[otheridx]
      end
      if diff == 0
        new.push @arr[selfidx]
        selfidx += 1
        otheridx += 1
      elsif diff < 0
        selfidx += 1
      elsif diff >0
        otheridx += 1
      end
    end

    return new
  end

  def -(other)
    new = []
    self.sort!
    other.sort!
    selfidx = 0
    otheridx = 0
    selflen = self.length
    otherlen = other.length

    while selfidx < selflen and otheridx < otherlen
      diff = 0
      if not @arr[selfidx].overlapping?(other[otheridx])
        diff = @arr[selfidx] <=> other[otheridx]
      end
      if diff == 0
        selfidx += 1
        otheridx += 1
      elsif diff < 0
        new.push @arr[selfidx]
        selfidx += 1
      elsif diff > 0
        otheridx += 1
      end
    end
    new += @arr[selfidx..-1]
  end

  def to_s
    @arr.join("\n")
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #$0 blastoutput condenseqoutput"
  STDERR.puts "format of blastout:"
  STDERR.puts "qseqid sseqid pident length qstart qend " +
    "sstart send evalue bitscore"
  exit 1
end

blastname = ARGV[0]
condenseqname = ARGV[1]

STDERR.puts "## sorting blasthits"
blastres = HitCollection.new(blastname).sort!
File.open('blastsorted', 'w') do |file|
  file.puts blastres
end

STDERR.puts "## sorting searchhits"
condres = HitCollection.new(condenseqname).sort!
File.open('searchsorted', 'w') do |file|
  file.puts condres
end

STDERR.puts "## &"
tp = condres & blastres
puts "## TP: #{tp.length}"
tp.each do |hit|
  puts hit
end

STDERR.puts "## search - blast"
fp = condres - blastres
puts "## FP: #{fp.length}"
fp.each do |hit|
  puts hit
end

STDERR.puts "## blast - search"
fn = blastres - condres
puts "## FN: #{fn.length}"
fn.each do |hit|
  puts hit
end
