module Fasta

require 'zlib'
# require 'codon2aa.rb'

class SequenceEntry
  def initialize(headline)
    h = headline.match(/^>(.*)\n/)
    if h
      @header = h[1]
    else
      STDERR.puts "#{$0}: illegal header #{header}"
      exit 1
    end
    @comment = Array.new()
    @seqoflines = Array.new()
  end
  def add_comment_line(line)
    @comment.push(line)
  end
  def add_sequence_line(line)
    @seqoflines.push(line.chomp.gsub(/\s/,""))
  end
  def set_header(s)
    @header = s
  end
  def write(rc, io, seq_line_width = 70)
    io.print ">#{@header}"
    @comment.each do |c|
      io.print c
    end
    s = @seqoflines.join("")
    if rc
      s = s.reverse
      s.tr!("ACGTacgt","TGCAtgca")
    end
    b = 0
    while b < s.size
      io.puts(s[b..b + seq_line_width - 1])
      b += seq_line_width
    end
  end
  def get_seqlength()
    lensum = 0
    @seqoflines.each do |line|
      lensum += line.length
    end
    return lensum
  end
  def get_sequence()
    return @seqoflines.join("")
  end
  def get_header()
    return @header
  end
end

=begin
  def write_open_reading_frame( offset, io, reverse=false )
    s = @seqoflines.join("")
    if reverse
      s = s.reverse
    end
    aaSequence = ""
    position = offset
    while position+3 < s.length    # enough bases for a codon
      begin
        codon = s[ position, 3 ]
        aa = codon2aa_11( codon )
        aaSequence += aa
      rescue => text
        io.puts text
        return
      end
      position += 3
    end
    puts aaSequence
  end
=end

def Fasta.read_multi_file(fname)  # function for module Fasta
  begin
    if fname.match(/\.gz$/)
      infp = Zlib::GzipReader.open(fname)
    else
      infp = File.open(fname,"r")
    end
  rescue => err
    STDERR.print "Could not open file \"#{fname}\": #{err}\n"
    exit 1
  end
  curr_entry = nil
  infp.each_line do |line|
    if line.match(/^>/)
      if not curr_entry.nil?
        yield curr_entry  # deliver current entry to iterator
      end
      curr_entry = SequenceEntry.new(line)  # create new sequence entry
    elsif line.match(/^;/)  # comment line
      curr_entry.add_comment_line(line)
    else
      curr_entry.add_sequence_line(line)
    end
  end
  if curr_entry.nil?
    STDERR.puts "#{$0}: assertion in f. #{__FILE__}, l. #{__LINE__} failed"
  end
  yield curr_entry
  infp.close
end

end # module Fasta
