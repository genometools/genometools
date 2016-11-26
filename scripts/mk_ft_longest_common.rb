#!/usr/bin/env ruby

def gen_access_raw(mode,structvar,pre)
  if mode == "bytes"
    return "*#{pre}ptr"
  elsif mode == "twobit"
    return "gt_twobitencoding_char_at_pos(\n" +
           " " * 30 + "#{structvar}->twobitencoding,\n" +
           " " * 30 + "#{pre}ptr)"
  elsif mode == "encseq"
    return "gt_encseq_get_encoded_char(#{structvar}->encseq,\n" +
           " " * 20 + "#{pre}ptr,\n" +
           " " * 20 + "GT_READMODE_FORWARD)"
  else
    return "gt_sequenceobject_esr_get(#{structvar},#{pre}ptr)"
  end
end

def gen_wildcard(wildcard)
  if wildcard
    return "cu == WILDCARD ||\n"
  else
    return ""
  end
end

def gen_compare(a_mode,b_mode,wildcard,complement)
  return "const GtUchar cu = #{gen_access_raw(a_mode,"useq","u")};\n" +
         " " * 8 + "if (" + gen_wildcard(wildcard) +
         (if wildcard then (" " * 12) else "" end) +
         "cu != " +
         (if complement then "GT_COMPLEMENTBASE(" else "" end) +
         "#{gen_access_raw(b_mode,"vseq","v")}" +
         (if complement then ")" else "" end) + ")\n" +
         " " * 10 + "break"
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

def gen_ptr_assign(mode,pre,left2right)
  if mode == "bytes"
    if left2right
      return "#{pre}seq->bytesequenceptr + #{pre}seq->offset + #{pre}start; " +
             "#{pre}step = 1"
    else
      return "#{pre}seq->bytesequenceptr + #{pre}seq->offset - #{pre}start; " +
             "#{pre}step = -1"
    end
  else
    if left2right
      if mode != "encseq_reader"
        return "#{pre}seq->offset + #{pre}start; #{pre}step = 1"
      else
        return "#{pre}start"
      end
    else
      if mode != "encseq_reader"
        return "#{pre}seq->offset - #{pre}start; #{pre}step = -1"
      else
        return "#{pre}start"
      end
    end
  end
end

def gen_ptr_set(mode,pre)
  if mode == "encseq_reader"
    return "#{pre}ptr = #{pre}start;"
  else
    return ["if (#{pre}seq->read_seq_left2right)",
            "{",
            "  #{pre}ptr = #{gen_ptr_assign(mode,pre,true)};",
            "} else",
            "{",
            "  #{pre}ptr = #{gen_ptr_assign(mode,pre,false)};",
            "}"].join("\n    ")
  end
end

def gen_var_decl(mode,pre)
  if mode == "bytes"
    return "const GtUchar *#{pre}ptr"
  else
    return "GtUword #{pre}ptr"
  end
end

def gen_step_decl(mode,pre)
  if mode != "encseq_reader"
    return "int #{pre}step;"
  else
    return ""
  end
end

def gen_ptr_incr(mode,pre)
  if mode == "encseq_reader"
    return "#{pre}ptr++"
  else
    return "#{pre}ptr += #{pre}step"
  end
end

def longestcommonfunc(a_mode,b_mode,wildcard)
 puts <<EOF
static GtUword #{gen_func_name(a_mode,b_mode,wildcard)}(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    #{gen_var_decl(a_mode,"u")};#{gen_step_decl(a_mode,"u")}
    #{gen_var_decl(b_mode,"v")};#{gen_step_decl(b_mode,"v")}
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    #{gen_ptr_set(a_mode,"u")}
    #{gen_ptr_set(b_mode,"v")}
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           #{gen_ptr_incr(a_mode,"u")}, #{gen_ptr_incr(b_mode,"v")}, matchlength++)
      {
        #{gen_compare(a_mode,b_mode,wildcard,true)};
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           #{gen_ptr_incr(a_mode,"u")}, #{gen_ptr_incr(b_mode,"v")}, matchlength++)
      {
        #{gen_compare(a_mode,b_mode,wildcard,false)};
      }
    }
    return matchlength;
  }
  return 0;
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
firstfunc = true
puts "\nGtLongestCommonFunc ft_longest_common_func_tab[] =\n{"
func_list.each_with_index do |func_name,idx|
  if firstwildcard.nil? and func_name.match(/_wildcard/)
    firstwildcard = idx
  end
  if firstfunc
    firstfunc = false
  else
    puts ","
  end
  print "  /* #{idx} */ #{func_name}"
end
puts "\n};"
puts "const int ft_longest_common_num_modes = #{modes.length};"
puts "const int ft_longest_common_func_first_wildcard = #{firstwildcard};"
