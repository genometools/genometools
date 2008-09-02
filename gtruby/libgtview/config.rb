#
# Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
require 'libgtcore/error'
require 'libgtview/color'

module GT
  extend DL::Importable
  gtdlload "libgt"

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

  extern "Style* style_new(bool, Error*)"
  extern "int style_load_file(Style*, Str*, Error*)"
  extern "int style_load_str(Style*, Str*, Error*)"
  extern "int style_to_str(const Style*, Str*, Error*)"
  extern "bool style_get_color(Style*, const char*, const char*, Color*)"
  extern "void style_set_color(Style*, const char*, const char*, Color*)"
  extern "bool style_get_str(const Style*, const char*, " +
                             "const char*, Str*)"
  extern "void style_set_str(Style*, const char*, const char*, Str*)"
  extern "bool style_get_num(const Style*, const char*, " +
                             "const char*, double*)"
  extern "void style_set_num(Style*, const char*, const char*, double)"
  extern "bool style_get_bool(const Style*, const char*, " +
                               "const char*, bool*)"
  extern "void style_set_bool(Style*, const char*, const char*, bool)"
  extern "void style_unset(Style*, const char*, const char*)"
  extern "void style_delete(Style*)"

  class Config
    attr_reader :config
    def initialize
      err = GT::Error.new()
      @config = GT.style_new(false, err.to_ptr)
      if not @config then GT.gterror(err) end
      @config.free = GT::symbol("style_delete", "0P")
    end

    def load_file(filename)
      err = GT::Error.new()
      str = GT::Str.new(filename)
      rval = GT.style_load_file(@config, str.to_ptr, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def load_str(str)
      err = GT::Error.new()
      str = GT::Str.new(str)
      rval = GT.style_load_str(@config, str.to_ptr, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def to_str()
      err = GT::Error.new()
      str = GT::Str.new(nil)
      if GT.style_to_str(@config, str.to_ptr, err.to_ptr) == 0
        str.to_s
      else
        GT.gterror(err)
      end
    end

    def clone
      cfg = GT::Config.new()
      str = self.to_str()
      cfg.load_str(str)
      cfg
    end

    def get_color(section, key)
      color = GT::Color.malloc
      if GT.style_get_color(@config, section, key, color)
        color
      else
        nil
      end
    end

    def set_color(section, key, color)
      GT.style_set_color(@config, section, key, color)
    end

    def get_cstr(section, key)
      str = GT::Str.new(nil)
      if GT.style_get_str(@config, section, key, str)
        str.to_s
      else
        nil
      end
    end

    def set_cstr(section, key, value)
      str = GT::Str.new(value)
      GT.style_set_str(@config, section, key, str)
    end

    def get_num(section, key)
      double = DoubleArg.malloc
      if GT.style_get_num(@config, section, key, double)
        double.val
      else
        nil
      end
    end

    def set_num(section, key, number)
      num = number.to_f
      GT.style_set_num(@config, section, key, num)
    end

    def get_bool(section, key)
      bool = BoolArg.malloc
      if GT.style_get_bool(@config, section, key, bool)
        bool.val
      else
        nil
      end
    end

    def set_bool(section, key, val)
      GT.style_set_bool(@config, section, key, val)
    end

    def unset(section, key)
      GT.style_unset(@config, section, key)
    end
  end
end
