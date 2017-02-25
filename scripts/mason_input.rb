require "scripts/fasta.rb"
require "scripts/polishing.rb"

def evaluate_cigarstring(cigarstring)
  costs = 0
  querylen = 0
  dblen = 0
  cigarstring.cigar_each do |multiplier,op|
    if op == 'D'
      dblen += multiplier
      costs += multiplier
    elsif op == 'I'
      querylen += multiplier
      costs += multiplier
    elsif op == 'X'
      dblen += multiplier
      querylen += multiplier
      costs += multiplier
    elsif op == 'M'
      dblen += multiplier
      querylen += multiplier
    else
      STDERR.puts "#{0}: unknown operator #{op} in cigar string"
      exit 1
    end
  end
  return costs, dblen, querylen
end

SampledSeqKeyvalues = Struct.new("SampledSeqKeyvalues",
                                 :costs, 
                                 :dblen, 
                                 :querylen, 
                                 :prefix_positive, 
                                 :suffix_positive, 
                                 :begin_pos,
                                 :strand,
                                 :identity)

def analyze_mason_header(header,polishing)
  dblen2 = nil
  sskv = SampledSeqKeyvalues.new(nil,nil,nil,nil,nil,nil,nil,nil)
  a = header.split(/\s/)
  1.upto(a.length-1).each do |idx|
    m = a[idx].match(/([A-Z_]*)=(.*)/)
    if not m
      STDERR.puts "#{$0}: cannot analyze #{header}"
      exit 1
    end
    key = m[1]
    value = m[2]
    if key == "CIGAR"
      sskv.costs, dblen2, sskv.querylen = evaluate_cigarstring(value)
      sskv.prefix_positive = polishing.prefix_positive?(value)
      sskv.suffix_positive = polishing.suffix_positive?(value)
    elsif key == "SAMPLE_SEQUENCE"
      sskv.dblen = value.length
    elsif key == "BEGIN_POS"
      sskv.begin_pos = value.to_i
    elsif key == "STRAND"
      sskv.strand = value
    end
  end
  if not sskv.costs.nil?
    if sskv.costs == 0
      sskv.identity = 100.0
    else
      sskv.identity = 100.0 - 
                      200.0 * sskv.costs.to_f/(sskv.dblen + sskv.querylen).to_f
    end
  end
  return sskv
end

def enumerate_valid_samples(readfile,minidentity)
  history = 60
  if minidentity.nil?
    error_rate = 1
  else
    error_rate = 100 - minidentity
  end
  polishing = Polishing.new(error_rate,history)
  Fasta.read_multi_file(readfile) do |seqentry|
    header = seqentry.get_header()
    sskv = analyze_mason_header(header,polishing)
    if sskv.prefix_positive and sskv.suffix_positive and 
       (minidentity.nil? or sskv.identity >= minidentity.to_f)
      newheader = "BEGIN_POS=#{sskv.begin_pos} " +
                  "STRAND=#{sskv.strand} " +
                  sprintf("IDENTITY=%.2f",sskv.identity)
      seqentry.set_header(newheader)
      yield seqentry
    end
  end
end
