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
    if self.overlapping?(other)
      return 0
    end
    if (@qseqid < other.qseqid)
      return -1
    elsif (@qseqid > other.qseqid)
      return 1
    end
    if (@sseqid < other.sseqid)
      return -1
    elsif (@sseqid > other.sseqid)
      return 1
    end
    if (@qstart < other.qstart)
      return -1
    elsif (@qstart > other.qstart)
      return 1
    end
    if (@sstart < other.sstart)
      return -1
    elsif (@sstart > other.sstart)
      return 1
    end
    if (@qend < other.qend)
      return -1
    elsif (@qend > other.qend)
      return 1
    end
    if (@send < other.send)
      return -1
    elsif (@send > other.send)
      return 1
    end
    return 0
  end

  def to_s
    "#@qseqid #@sseqid #@pident #@length #@qstart #@qend #@sstart #@ssend" +
      " #@evalue #@bitscore"
  end

  def overlapping?(other)
    if (@qseqid == other.qseqid and @sseqid == other.sseqid)
      overlapstart = [@qstart, other.qstart].max
      overlapend = [@qend, other.qend].min
      overlaplength = overlapend - overlapstart
      if ((@length * 0.8) <= overlaplength and
        (other.length * 0.8) <= overlaplength)
        return true
      end
    end
    return false
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
      case @arr[selfidx] <=> other[otheridx]
      when 0
        new.push @arr[selfidx]
        selfidx += 1
        otheridx += 1
      when -1
        selfidx += 1
      when 1
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
      case @arr[selfidx] <=> other[otheridx]
      when 0
        selfidx += 1
        otheridx += 1
      when -1
        new.push @arr[selfidx]
        selfidx += 1
      when 1
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
  STDERR.puts "Usage: #$0 blastoutput condenseroutput"
  STDERR.puts "format of blastout:"
  STDERR.puts "qseqid sseqid pident length qstart qend " +
    "sstart send evalue bitscore"
  exit 1
end

blastname = ARGV[0]
condensername = ARGV[1]

STDERR.puts "## sorting blasthits"
blastres = HitCollection.new(blastname).sort!
File.open('blastsorted', 'w') do |file|
  file.puts blastres
end

STDERR.puts "## sorting searchhits"
condres = HitCollection.new(condensername).sort!
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
