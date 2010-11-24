#!/usr/bin/env ruby
#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

class LTRAnnotation
  attr_accessor :indexm, :indexs, :data, :file, :commentlines

  def initialize
    @data = []
    @commentlines = []
    @indexm = {}
    @indexs = {}
    @file = ""
  end

  # The 'load_from_*' functions require at least that the first two items on
  # a line represent a LTR retrotransposon start and end. How the rest of each
  # line is read and processed depends on the file type, according to each
  # load function.

  def load_from_ltrharvest(file)
    File.open(file) do |ltr|
    @file = file
      ltr.each_line do |line|
        if line[0].chr != "#" then
          (sltr,eltr,lltr,sleft,eleft,lleft,sright,eright,lright) = \
            line.split(" ").collect{|i| i.strip}
          smiddle = (eleft.to_i+1).to_s
          newitem = {:sltr => sltr.to_i,
                    :eltr => eltr.to_i,
                    :lltr => lltr.to_i,
                    :sleft => sleft.to_i,
                    :eleft => eleft.to_i,
                    :lleft => lleft.to_i,
                    :sright => sright.to_i,
                    :eright => eright.to_i,
                    :lright => lright.to_i}
          @data.push(newitem)
          @indexm[smiddle] = newitem
          @indexs[sltr] = newitem
        else
          @commentlines.push(line.chop)
        end
      end
    end
  end

  def load_from_berkeley(file)
    File.open(file) do |anno|
    @file=file
      anno.each_line do |line|
        if line[0].chr != "#" then
          (sltr, eltr, fbti_id, family_name, fbgn_id, chromosome_arm, release)\
             = line.split(" ").collect{|i| i.strip}
          newitem = {:sltr => sltr.to_i,
                    :eltr => eltr.to_i,
                    :length => eltr.to_i - sltr.to_i + 1,
                    :fbti_id => fbti_id,
                    :family_name => family_name,
                    :fbgn_id => fbgn_id,
                    :chromosome_arm => chromosome_arm,
                    :release => release}
          @data.push(newitem)
          @indexs[sltr] = newitem
        else
          @commentlines.push(line.chop)
        end
      end
    end
  end

  def load_from_custom(file)
    File.open(file) do |anno|
    @file=file
      anno.each_line do |line|
        if line[0].chr != "#" then
          (sltr, eltr, lltr, lname, ltype, lstrand)\
             = line.split(" ").collect{|i| i.strip}
          newitem = {:sltr => sltr.to_i,
                    :eltr => eltr.to_i,
                    :lltr => lltr.to_i,
                    :lname => lname,
                    :ltype => ltype,
                    :lstrand => lstrand}
          @data.push(newitem)
          @indexs[sltr] = newitem
        else
          @commentlines.push(line.chop)
        end
      end
    end
  end

  def print_content
    @data.each do |item|
      puts item.inspect
    end
  end

  def merge_hits(ref_hit, pred_hit)
    return {:s_pred => pred_hit[:sltr],
            :e_pred => pred_hit[:eltr],
            :s_ref => ref_hit[:sltr],
            :e_ref => ref_hit[:eltr],
            :length => ref_hit[:length],
            :fbti_id => ref_hit[:fbti_id],
            :family_name => ref_hit[:family_name],
            :fbgn_id => ref_hit[:fbgn_id],
            :chromosome_arm => ref_hit[:chromosome_arm],
            :release => ref_hit[:release]
            }
  end

  #compare a LTRharvest annotation with a reference annotation
  def compare(refAnno, difference=0, shift=0, without_ltrs=false)
    tp_list           = []
    htp_startpos_list = []
    htp_endpos_list   = []
    fp_list           = []
    fn_list           = []
    foundlist         = []
    nof_fn            = 0
    @data.each do |item|
      found = false
      item[:sltr] += shift
      item[:eltr] += shift
      if without_ltrs then
        item[:sltr] += item[:lleft]
        item[:eltr] -= item[:lright]
      end
      refAnno.data.each do |ref_item|
        p_start = false
        p_end = false
        if (ref_item[:sltr] - difference <= item[:sltr] && \
            item[:sltr] <= ref_item[:sltr] + difference) then
          p_start = true
        end
        if (ref_item[:eltr] - difference <= item[:eltr] && \
            item[:eltr] <= ref_item[:eltr] + difference) then
          p_end = true
        end
        if (p_start and not p_end) then
          htp_startpos_list.push(merge_hits(ref_item, item))
          foundlist.push([ref_item[:sltr],ref_item[:eltr]])
          found = true
        elsif (p_end and not p_start) then
          htp_endpos_list.push(merge_hits(ref_item, item))
          foundlist.push([ref_item[:sltr],ref_item[:eltr]])
          found = true
        elsif (p_start and p_end) then
          tp_list.push(merge_hits(ref_item, item))
          foundlist.push([ref_item[:sltr],ref_item[:eltr]])
          found = true
        end
      end
      if !found then
        fp_list.push(item)
      end
    end
    fn_list = refAnno.data.reject do |ref_item|
      foundlist.include?([ref_item[:sltr],ref_item[:eltr]])
    end
    yield tp_list, htp_startpos_list, htp_endpos_list,\
          fp_list, fn_list, difference
  end

  #output comparison results in old evalscript format
  def compare_output(refAnno, difference)
    compare(refAnno, 20) do |tp_list, htp_startpos_list, htp_endpos_list,\
                           fp_list, fn_list, difference|
      puts "# compared start- and endpositions from \n#  file \"#{@file}\" and\n#"\
           + "  file \"#{refAnno.file}\"\n#"
      puts "# parameters from LTRharvest (comments from\n#   file \"#{@file}\"):"
      @commentlines.each do |line|
        puts line
      end
      puts "#\n# read #{refAnno.data.length} annotations of" +\
         " LTR retrotransposons from\n" +\
         "#   file #{refAnno.file} \n#"
      puts("# results: each row = 1 LTR-pair\n" +\
         "#          startpos = 5'-boundary (of leftLTR)\n" +\
         "#          endpos   = 3'-boundary (of rightLTR)\n#\n" +\
         "# start- and endpositions may vary +/-#{difference}" +\
         "bp.\n#")

      puts "# true  positive(s) (TP) = #{tp_list.length}"
      puts "# half true positive(s) (TP) = #{htp_startpos_list.length\
                                           + htp_endpos_list.length}"
      puts "# false positive(s) (FP) = #{fp_list.length}"
      puts "# false negative(s) (FN) = #{fn_list.length}"
      puts "#"
#      puts "# duplicated TP predictions = "  #?
      puts "#"
      #print TPs
      puts "# #{tp_list.length} true positive(s) (TP)\n" +\
           "# startposition  endposition  " +\
           "+/-startposition_berkeley  +/-endposition_berkeley  " +\
           "length  FBti_id  family_name  FBgn_id  chromosome_arm  release"
      tp_list.sort{|i1, i2| i1[:s_pred] <=> i2[:s_pred]}.each do |tp|
        puts "#{tp[:s_pred]} #{tp[:e_pred]} #{tp[:s_ref]-tp[:s_pred]} "\
             +"#{tp[:e_ref]-tp[:e_pred]} #{tp[:length]} #{tp[:fbti_id]} "\
             +"#{tp[:family_name]} #{tp[:fbgn_id]} #{tp[:chromosome_arm]} "\
             +"#{tp[:release]}"
      end
      #print hTPs
      puts "# #{htp_startpos_list.length} half true positive(s) (only first position correct)\n" +\
           "# startposition  false_endposition  "+\
           "+/-startposition_berkeley  berkeley_endposition  "+\
           "length  FBti_id  family_name  FBgn_id  chromosome_arm  release"
      htp_startpos_list.sort{|i1, i2| i1[:s_pred] <=> i2[:s_pred]}.each do |htp|
        puts "#{htp[:s_pred]} #{htp[:e_pred]} #{htp[:s_ref]-htp[:s_pred]} "\
             +"#{htp[:e_ref]} #{htp[:length]} #{htp[:fbti_id]} "\
             +"#{htp[:family_name]} #{htp[:fbgn_id]} #{htp[:chromosome_arm]} "\
             +"#{htp[:release]}"
      end
      puts "# #{htp_endpos_list.length} half true positive(s) (only second position correct)\n" +\
            "# false_startposition  endposition  "+\
            "berkeley_startposition  +/-endposition_berkeley"+\
            "length  FBti_id  family_name  FBgn_id  chromosome_arm  release"
      htp_endpos_list.sort{|i1, i2| i1[:s_pred] <=> i2[:s_pred]}.each do |htp|
        puts "#{htp[:s_pred]} #{htp[:e_pred]} #{htp[:s_ref]} "\
             +"#{htp[:e_ref]-htp[:e_pred]} #{htp[:length]} #{htp[:fbti_id]} "\
             +"#{htp[:family_name]} #{htp[:fbgn_id]} #{htp[:chromosome_arm]} "\
             +"#{htp[:release]}"
      end
      #print FPs
      puts "# #{fp_list.length} false positive(s) (FP)\n" +\
           "# false_startposition  false_endposition  length  " +\
           "startleftLTR  endleftLTR  lengthleftLTR  " +\
          "startrightLTR  endrightLTR  lengthrightLTR"
      fp_list.sort{|i1, i2| i1[:sltr] <=> i2[:sltr]}.each do |fp|
        puts "#{fp[:sltr]} #{fp[:eltr]} #{fp[:lltr]} "\
             +"#{fp[:sleft]} #{fp[:eleft]} #{fp[:lleft]} "\
             +"#{fp[:sright]} #{fp[:eright]} #{fp[:lright]}"
      end
      #print FNs
      puts "# #{fn_list.length} false negative(s) (FN)\n" +\
           "# not predicted startposition, not predicted endposition"
      fn_list.sort{|i1, i2| i1[:sltr] <=> i2[:sltr]}.each do |fn|
        puts "#{fn[:sltr]} #{fn[:eltr]} #{fn[:length]} #{fn[:fbti_id]} "\
             +"#{fn[:family_name]} #{fn[:fbgn_id]} #{fn[:chromosome_arm]} "\
             +"#{fn[:release]}"
      end
    end
  end
end

#mixin for arbitrary decimal rounding
class Float
  def round_to(x)
    (self * 10**x).round.to_f / 10**x
  end

  def ceil_to(x)
    (self * 10**x).ceil.to_f / 10**x
  end

  def floor_to(x)
    (self * 10**x).floor.to_f / 10**x
  end
end
