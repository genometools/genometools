#
# Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c)      2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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

require 'gtdlload'
require 'gthelper'
require 'core/error'
require 'annotationsketch/color'

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  # looks weird, but apparently the way to pass proper pointers to external
  # functions
  DoubleArg = struct [
    "double val"
  ]

  typealias "bool", "ibool"

  # same as above
  BoolArg = struct [
    "bool val"
  ]

  # a NULL pointer
  NULL = DL::PtrData.new(0)

  extern "GT_Style* gt_style_new(bool, GT_Error*)"
  extern "int gt_style_load_file(GT_Style*, const char*, GT_Error*)"
  extern "int gt_style_load_str(GT_Style*, GT_Str*, GT_Error*)"
  extern "int gt_style_to_str(const GT_Style*, GT_Str*, GT_Error*)"
  extern "bool gt_style_get_color(GT_Style*, const char*, const char*, GT_Color*, " +
                                 "GenomeNode*)"
  extern "void gt_style_set_color(GT_Style*, const char*, const char*, GT_Color*)"
  extern "bool gt_style_get_str(const GT_Style*, const char*, " +
                               "const char*, GT_Str*, GT_GenomeNode*)"
  extern "void gt_style_set_str(GT_Style*, const char*, const char*, GT_Str*)"
  extern "bool gt_style_get_num(const GT_Style*, const char*, " +
                               "const char*, double*, GT_GenomeNode*)"
  extern "void gt_style_set_num(GT_Style*, const char*, const char*, double)"
  extern "bool gt_style_get_bool(const GT_Style*, const char*, " +
                                "const char*, bool*, GT_GenomeNode*)"
  extern "void gt_style_set_bool(GT_Style*, const char*, const char*, bool)"
  extern "void gt_style_unset(GT_Style*, const char*, const char*)"
  extern "void gt_style_delete(GT_Style*)"

  class Style
    attr_reader :style

    def initialize
      err = GT::Error.new()
      @style = GT.gt_style_new(false, err.to_ptr)
      if not @style then GT.gterror(err) end
      @style.free = GT::symbol("gt_style_delete", "0P")
    end

    def load_file(filename)
      err = GT::Error.new()
      rval = GT.gt_style_load_file(@style, filename, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def load_str(str)
      err = GT::Error.new()
      str = GT::Str.new(str)
      rval = GT.gt_style_load_str(@style, str.to_ptr, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def to_str()
      err = GT::Error.new()
      str = GT::Str.new(nil)
      if GT.gt_style_to_str(@style, str.to_ptr, err.to_ptr) == 0
        str.to_s
      else
        GT.gterror(err)
      end
    end

    def clone
      sty = GT::Style.new()
      str = self.to_str()
      sty.load_str(str)
      sty
    end

    def get_color(section, key, gn = GT::NULL)
      color = GT::Color.malloc
      if GT.gt_style_get_color(@style, section, key, color, gn)
        color
      else
        nil
      end
    end

    def set_color(section, key, color)
      GT.gt_style_set_color(@style, section, key, color)
    end

    def get_cstr(section, key, gn = GT::NULL)
      str = GT::Str.new(nil)
      if GT.gt_style_get_str(@style, section, key, str, gn)
        str.to_s
      else
        nil
      end
    end

    def set_cstr(section, key, value)
      str = GT::Str.new(value)
      GT.gt_style_set_str(@style, section, key, str)
    end

    def get_num(section, key, gn = GT::NULL)
      double = DoubleArg.malloc
      if GT.gt_style_get_num(@style, section, key, double, gn)
        double.val
      else
        nil
      end
    end

    def set_num(section, key, number)
      num = number.to_f
      GT.gt_style_set_num(@style, section, key, num)
    end

    def get_bool(section, key, gn = GT::NULL)
      bool = BoolArg.malloc
      if GT.gt_style_get_bool(@style, section, key, bool, gn)
        bool.val
      else
        nil
      end
    end

    def set_bool(section, key, val)
      GT.gt_style_set_bool(@style, section, key, val)
    end

    def unset(section, key)
      GT.gt_style_unset(@style, section, key)
    end
  end
end
