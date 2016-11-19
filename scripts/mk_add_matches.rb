#!/usr/bin/env ruby

require "erb"

class ERBContext
  def initialize(hash)
    hash.each_pair do |key, value|
      instance_variable_set('@' + key.to_s, value)
    end
  end
  def get_binding
    binding
  end
end

def gen_access_expr(posvar,structvar)
  return "#{structvar}->read_seq_left2right\n" +
         " " * 45 + "? #{structvar}->offset + #{posvar}\n" +
         " " * 45 + ": #{structvar}->offset - #{posvar}"
end

def gen_access_raw(mode,posvar,structvar)
  if mode == "bytes"
    return "#{structvar}->bytesequenceptr[\n" +
           " " * 20 + gen_access_expr(posvar,structvar) + "\n" +
           " " * 20 + "]"
  elsif mode == "twobit"
    return "gt_twobitencoding_char_at_pos(\n" +
           " " * 20 + "#{structvar}->twobitencoding,\n" +
           " " * 20 + gen_access_expr(posvar,structvar) + "\n" +
           " " * 20 + ")"
  elsif mode == "encseq"
    return "gt_encseq_get_encoded_char(#{structvar}->encseq,\n" +
           " " * 20 + gen_access_expr(posvar,structvar) + ",\n" +
           " " * 20 + "GT_READMODE_FORWARD)"
  end
end

def gen_access(mode,isleft)
  if isleft
    return gen_access_raw(mode,"upos","useq")
  else
    s = gen_access_raw(mode,"vpos","vseq")
    return "(vseq->dir_is_complement\n" +
           " " * 10 + "? GT_COMPLEMENTBASE(#{s})\n" +
           " " * 10 + ": #{s})"
  end
end

def gen_special(special)
  if special 
    return "!ISSPECIAL(cu) || " 
  else 
    return ""
  end
end

def gen_suffix(special)
  if special
    return "_special"
  else
    return ""
  end
end

def cut_out(s)
  width = 80
  line = String.new()
  puts "code=#{s}"
  s.scan(/(\S+)/) do |item|
    puts "found #{item}"
    if line.length + item[1].length <= width
      line = line + item[1]
    else
      puts line
      line.clear
    end
  end
  if not line.empty?
    puts line
  end
end

def longestcommonfunc(a_mode,b_mode,special)
 erb_context = ERBContext.new({:a_mode => a_mode,
                               :b_mode => b_mode,
                               :special => special})
 header = ERB.new <<'EOF'
GtUword ft_longest_common_<%= @a_mode%>_<%= @b_mode%><%= gen_suffix(@special)%>(const GtFtFrontvalue *midfront,
                                      const GtFtFrontvalue *fv,
                                      GtFtSequenceObject *useq,
                                      GtFtSequenceObject *vseq)
{
  GtUword upos, vpos;

  for (upos = fv->row, vpos = fv->row + GT_FRONT_DIAGONAL(fv);
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = <%= gen_access(@a_mode,true)%>;

    if (<%= gen_special(@special)%>cu != <%= gen_access(@b_mode,false)%>)
    {
      break;
    }
  }
  return upos - fv->row;
}
EOF
  puts header.result(erb_context.get_binding)
end

first = true
modes = ["bytes","twobit","encseq"]
[true,false].each do |special|
  modes.each do |a_mode|
    modes.each do |b_mode|
     if first
       first = false
     else
       puts ""
     end
     longestcommonfunc(a_mode,b_mode,special)
    end
  end
end
