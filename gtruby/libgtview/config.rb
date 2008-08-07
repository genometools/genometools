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

  extern "Config* config_new(bool, Error*)"
  extern "int config_load_file(Config*, Str*, Error*)"
  extern "void config_get_color(Config*, const char*, const char*, Color*)"
  extern "void config_set_color(Config*, const char*, const char*, Color*)"
  extern "bool config_get_str(const Config*, const char*, " +
                             "const char*, Str*)"
  extern "void config_set_str(const Config*, const char*, " +
                             "const char*, Str*)"
  extern "bool config_get_num(const Config*, const char*, " +
                             "const char*, double*)"
  extern "void config_set_num(const Config*, const char*, " +
                             "const char*, double)"
  extern "bool config_get_bool(const Config*, const char*, " +
                               "const char*, bool*)"
  extern "void config_set_bool(const Config*, const char*, " +
                              "const char*, bool)"
  extern "bool config_get_cstr_list(const Config*, const char*, " +
                                   "const char*, StrArray*)"
  extern "void config_set_cstr_list(const Config*, const char*, " +
                                   "const char*, StrArray*)"
  extern "void config_delete(Config*)"

  class Config
    attr_reader :config
    def initialize
      err = GT::Error.new()
      @config = GT.config_new(false, err.to_ptr)
      if not @config then GT.gterror(err) end
      @config.free = GT::symbol("config_delete", "0P")
    end

    def load_file(filename)
      err = GT::Error.new()
      str = GT::Str.new(filename)
      rval = GT.config_load_file(@config, str.to_ptr, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def get_color(section, key)
      color = GT::Color.malloc
      GT.config_get_color(@config, section, key, color)
      color
    end

    def set_color(section, key, color)
      GT.config_set_color(@config, section, key, color)
    end

    def get_cstr(section, key)
      str = GT::Str.new(nil)
      if GT.config_get_str(@config, section, key, str)
        str.to_s
      else
        "undefined"
      end
    end

    def set_cstr(section, key, value)
      str = GT::Str.new(value)
      GT.config_set_str(@config, section, key, str)
    end

    def get_num(section, key)
      double = DoubleArg.malloc
      if GT.config_get_num(@config, section, key, double)
        double.val
      else
        -9999.99
      end
    end

    def set_num(section, key, number)
      num = number.to_f
      GT.config_set_num(@config, section, key, num)
    end

    def get_bool(section, key)
      bool = BoolArg.malloc
      if GT.config_get_bool(@config, section, key, bool)
        bool.val
      else
        false
      end
    end

    def set_bool(section, key, val)
      GT.config_set_bool(@config, section, key, val)
    end

    def get_cstr_list(section, key)
      strarr = GT::StrArray.new()
      GT.config_get_cstr_list(@config, section, key, strarr.strarray)
      strarr.to_a
    end

    def set_cstr_list(section, key, list)
      strarray = GT::StrArray.new()
      strarray.add_list(list)
      GT.config_set_cstr_list(@config, section, key, strarray.strarray)
    end
  end
end
