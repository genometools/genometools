#!/usr/bin/env ruby

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
  else
    return "gt_sequenceobject_esr_get(#{structvar},#{posvar})"
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

def gen_wildcard(wildcard)
  if wildcard
    return "(cu == WILDCARD) || "
  else
    return ""
  end
end

def gen_suffix(wildcard)
  if wildcard
    return "_wildcard"
  else
    return ""
  end
end

def gen_func_name(a_mode,b_mode,wildcard)
  return "ft_longest_common_#{a_mode}_#{b_mode}#{gen_suffix(wildcard)}"
end

def longestcommonfunc(a_mode,b_mode,wildcard)
 puts <<EOF
static GtUword #{gen_func_name(a_mode,b_mode,wildcard)}(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = #{gen_access(a_mode,true)};

    if (#{gen_wildcard(wildcard)}cu != #{gen_access(b_mode,false)})
    {
      break;
    }
  }
  return upos - ustart;
}
EOF
end

first = true
func_list = Array.new()
modes = ["twobit","encseq_reader","encseq","bytes"]
[false,true].each do |wildcard|
  modes.each do |a_mode|
    modes.each do |b_mode|
     if first
       first = false
     else
       puts ""
     end
     longestcommonfunc(a_mode,b_mode,wildcard)
     func_list.push(gen_func_name(a_mode,b_mode,wildcard))
    end
  end
end

firstwildcard = nil
puts "\nGtLongestCommonFunc ft_longest_common_func_tab[] =\n{"
print "  /* 0 */ ft_longest_common_all"
func_list.each_with_index do |func_name,idx|
  if firstwildcard.nil? and func_name.match(/_wildcard/)
    firstwildcard = idx
  end
  print ",\n  /* #{idx+1} */ #{func_name}"
end
puts "\n};"
puts "const int ft_longest_common_num_modes = #{modes.length};"
puts "const int ft_longest_common_func_first_wildcard = #{firstwildcard};"
