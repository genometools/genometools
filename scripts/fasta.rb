module Fasta

class SequenceEntry
  def initialize(header)
    @header = header
    @comment = []
    @seqoflines = []
  end
  def add_comment_line(line)
    @comment.push(line)
  end
  def add_sequence_line(line)
    line.chomp!
    line.gsub!(/\s/,"")
    @seqoflines.push(line)
  end
  def write(rc, io, seq_line_width = 70)
    io.print @header
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
      lensum+=line.length
    end
    return lensum
  end
  def get_sequence()
    return @seqoflines.join("")
  end
  def get_header()
    return @header.chomp
  end
end

def Fasta.read_multi_file(fname)  # function for module Fasta
  begin
    fl = File.open(fname,"r") 
  rescue => err
    STDERR.print "Could not open file \"#{fname}\": #{err}\n"
    exit 1
  end
  curr_entry = nil
  fl.each_line do |l|
    if l[0] == ?>
      if curr_entry != nil 
         yield curr_entry  # deliver current entry to iterator
      end 
      curr_entry = SequenceEntry.new(l)  # create new sequence entry
    elsif l[0] == ?;   # comment line
      curr_entry.add_comment_line(l)
    else
      curr_entry.add_sequence_line(l)
    end
  end
  if curr_entry then 
    yield curr_entry 
  end  # for the last sequence entry
end

end # module Fasta
