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

def gen_compare(a_mode,b_mode,wildcard,complement)
  access_v = gen_access_raw(b_mode,"vseq","v")
  splitter = if access_v.length > 10 then ("\n" + " " * 12) else " " end
  cmp_expr = (if complement then ("GT_COMPLEMENTBASE(") else "" end) +
             "#{access_v}" +
             (if complement then ")" else "" end) + ")\n" +
             " " * 10 + "break"
  if wildcard
    return "const GtUchar cu = #{gen_access_raw(a_mode,"useq","u")};\n" +
           " " * 8 + "if (cu == WILDCARD ||" + splitter + "cu !=" + splitter +
           cmp_expr
  else
    return "if (#{gen_access_raw(a_mode,"useq","u")} !=" + splitter +
           cmp_expr
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
  return ["\n    if (#{pre}seq->read_seq_left2right)",
          "{",
          "  #{pre}ptr = #{gen_ptr_assign(mode,pre,true)};",
          "} else",
          "{",
          "  #{pre}ptr = #{gen_ptr_assign(mode,pre,false)};",
          "}"].join("\n    ")
end

def gen_minsub(pre)
  other = if pre == "u" then "v" else "u" end
  return ["\n    GtUword minsubstringlength = #{other}start + #{pre}seq->substringlength - #{pre}start;",
          "if (#{other}seq->substringlength < minsubstringlength)",
          "{",
          "  minsubstringlength = #{other}seq->substringlength;",
          "}"].join("\n    ")
end

def gen_minsubstringlength_decl(a_mode,b_mode)
  minmatch_decl = ["\n    GtUword minsubstringlength = useq->substringlength - ustart;",
                     "if (vseq->substringlength < useq->insubstringlength)",
                     "{",
                     "  minsubstringlength = vseq->substringlength - vstart;",
                     "}"].join("\n    ") +
            gen_ptr_set(a_mode,"u") +
            gen_ptr_set(b_mode,"v")
  if a_mode == "encseq_reader"
    if b_mode == "encseq_reader"
      return "GtUword uptr = ustart, vptr = vstart;" +
             gen_minsub("v")
    elsif b_mode == "bytes"
      return "GtUword uptr = ustart; const GtUchar *vptr; int vstep;" +
             gen_minsub("v") +
             gen_ptr_set(b_mode,"v")
    else
      return "GtUword uptr = ustart, vptr; int vstep;" +
             gen_minsub("v") +
             gen_ptr_set(b_mode,"v")
    end
  elsif b_mode == "encseq_reader"
    if a_mode == "bytes"
      return "const GtUchar *uptr; int ustep; GtUword vptr = vstart;" +
             gen_minsub("u") +
             gen_ptr_set(a_mode,"u")
    else
       return "GtUword uptr; int ustep; GtUword vptr = vstart;" +
              gen_minsub("u") +
              gen_ptr_set(a_mode,"u")
    end
  else
    minmatch_decl = ["    GtUword minsubstringlength = useq->substringlength - ustart,",
                     "matchlength = 0;",
                     "if (vseq->substringlength - vstart < minsubstringlength)",
                     "{",
                     "  minsubstringlength = vseq->substringlength - vstart;",
                     "}"].join("\n    ") +
            gen_ptr_set(a_mode,"u") +
            gen_ptr_set(b_mode,"v")

    if a_mode == "bytes"
      if b_mode == "bytes"
        return "const GtUchar *uptr, *vptr; int ustep, vstep;\n" +
               minmatch_decl
      else
        return "const GtUchar *uptr; GtUword vptr; int ustep, vstep;\n" +
               minmatch_decl
      end
    else
      if b_mode == "bytes"
        return "GtUword uptr; const GtUchar *vptr; int ustep, vstep;\n" +
               minmatch_decl
      else
        return "GtUword uptr, vptr; int ustep, vstep;\n" +
               minmatch_decl
      end
    end
  end
end

def gen_ptr_incr(mode,pre)
  if mode == "encseq_reader"
    return "#{pre}ptr++"
  else
    return "#{pre}ptr += #{pre}step"
  end
end

def gen_matchlength_inc(a_mode,b_mode)
  if a_mode == "encseq_reader" or b_mode == "encseq_reader"
    return ""
  else
    return "matchlength++;"
  end
end

def gen_smaller(a_mode,b_mode)
  if a_mode == "encseq_reader"
    return "uptr"
  elsif b_mode == "encseq_reader"
    return "vptr"
  else
    return "matchlength"
  end
end

def gen_return_matchlength(a_mode,b_mode)
  if a_mode == "encseq_reader"
    return "uptr - ustart"
  elsif b_mode == "encseq_reader"
    return "vptr - vstart"
  else
    return "matchlength"
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
    #{gen_minsubstringlength_decl(a_mode,b_mode)}
    if (vseq->dir_is_complement)
    {
      do
      {
        #{gen_compare(a_mode,b_mode,wildcard,true)};
        #{gen_ptr_incr(a_mode,"u")};
        #{gen_ptr_incr(b_mode,"v")};#{gen_matchlength_inc(a_mode,b_mode)}
      } while (#{gen_smaller(a_mode,b_mode)} < minsubstringlength);
    } else
    {
      do
      {
        #{gen_compare(a_mode,b_mode,wildcard,false)};
        #{gen_ptr_incr(a_mode,"u")};
        #{gen_ptr_incr(b_mode,"v")};#{gen_matchlength_inc(a_mode,b_mode)}
      } while (#{gen_smaller(a_mode,b_mode)} < minsubstringlength);
    }
    return #{gen_return_matchlength(a_mode,b_mode)};
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
