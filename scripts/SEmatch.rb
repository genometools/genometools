class SEmatch
  def initialize(matchfile)
    @matchfile = matchfile
  end
  def each()
    begin
      @fp = File.new(@matchfile)
    rescue => err
      STDERR.puts "#{$0}: cannot open #{matchfile}"
      exit 1
    end
    @options = nil
    @fields = nil
    @runtime = nil
    @spacepeak = nil
    @fp.each_line do |line|
      if line.match(/^#/)
        if m = line.match(/^# Options:(.*)$/)
          @options = m[1]
        elsif m = line.match(/^# Fields:\s(.*)$/)
          @fields = m[1].gsub(/\./,"_").split(/, /)
          @fields.map!{|f| f.gsub(/\s/,"_")}
        end
        if m = line.match(/^# TIME.*\s(\S+)$/)
          @runtime = m[1].to_f
        elsif m = line.match(/^# combined space peak in megabytes:\s(\S+)$/)
          @spacepeak = m[1].to_f
        end
      else
        if @fields.nil?
          STDERR.puts "#{$0}: fields: not defined"
          exit 1
        end
        value_hash = Hash.new()
        line.gsub(/^\s+/,"").split(/\s+/).each_with_index do |value,idx|
          if value.match(/^\d+\.\d+$/)
            value_hash[@fields[idx].to_sym] = value.to_f
          elsif value.match(/^\d+$/)
            value_hash[@fields[idx].to_sym] = value.to_i
          else
            value_hash[@fields[idx].to_sym] = value
          end
        end
        if value_hash.has_key?(:s_start)
          if value_hash.has_key?(:s_len)
            value_hash[:s_end] = value_hash[:s_start] + value_hash[:s_len] - 1
          elsif value_hash.has_key?(:s_end)
            value_hash[:s_len] = value_hash[:s_end] - value_hash[:s_start] + 1
          else
            STDERR.puts "#{$0}: either length of match on subject or end " +
                        "position must be given"
            exit 1
          end
        else
          STDERR.puts "#{$0}: start of match on subject must be given"
          exit 1
        end
        if value_hash.has_key?(:q_start)
          if value_hash.has_key?(:q_len)
            value_hash[:q_end] = value_hash[:q_start] + value_hash[:q_len] - 1
          elsif value_hash.has_key?(:q_end)
            value_hash[:q_len] = value_hash[:q_end] - value_hash[:q_start] + 1
          else
            value_hash[:q_len] = value_hash[:s_len]
          end
        else
          STDERR.puts "#{$0}: start of match on query must be given"
          exit 1
        end
        value_hash[:origline] = line.chomp
        yield value_hash
      end
    end
  end
  def spacepeak_get()
    return @spacepeak
  end
  def runtime_get()
    return @runtime
  end
end

def match_is_identical(m0,m1)
  if m0[:s_seqnum] != m1[:s_seqnum] or m0[:q_seqnum] != m1[:q_seqnum]
    STDERR.puts "expect same sequence numbers"
    exit 1
  end
  if [m0[:s_start],m0[:s_end],m0[:q_start],m0[:q_end]] ==
     [m1[:s_start],m1[:s_end],m1[:q_start],m1[:q_end]]
    return true
  end
  return false
end

def coords_contained_in(start0,end0,start1,end1)
  if start1 <= start0 and end0 <= end1
    return true
  end
  return false
end

def match_proper_contained_in(m0,m1)
  if m0[:s_seqnum] != m1[:s_seqnum] or m0[:q_seqnum] != m1[:q_seqnum]
    STDERR.puts "expect same sequence numbers"
    exit 1
  end
  if coords_contained_in(m0[:s_start],m0[:s_end],m1[:s_start],m1[:s_end]) and
     coords_contained_in(m0[:q_start],m0[:q_end],m1[:q_start],m1[:q_end]) and
     not match_is_identical(m0,m1)
    return true
  end
  return false
end

def coords_overlap_size(start0,end0,start1,end1)
  if start0 > start1
    STDERR.puts "start0=#{start0} > #{start1}=start1 not expected"
    exit 1
  end
  if end0 < start1
    return 0
  elsif end0 < end1
    return end0 - start1 + 1
  else
    return end1 - start1 + 1
  end
end

def matches_overlap(m0,m1)
  ovl = 0
  if m0[:s_start] <= m1[:s_start]
    ovl += coords_overlap_size(m0[:s_start],m0[:s_end],m1[:s_start],m1[:s_end])
  else
    ovl += coords_overlap_size(m1[:s_start],m1[:s_end],m0[:s_start],m0[:s_end])
  end
  if m0[:q_start] <= m1[:q_start]
    ovl += coords_overlap_size(m0[:q_start],m0[:q_end],m1[:q_start],m1[:q_end])
  else
    ovl += coords_overlap_size(m1[:q_start],m1[:q_end],m0[:q_start],m0[:q_end])
  end
  len = (m0[:s_len] + m0[:q_len])/2 + (m1[:s_len] + m1[:q_len])/2
  return (ovl.to_f/len.to_f).round
end
