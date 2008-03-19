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
  typealias "bool", "ibool"
  extern "Config* config_new(bool, Error*)"
  extern "int config_load_file(Config*, Str*, Error*)"
  extern "void config_get_colorptr(Config*, Color*, const char*)"
  extern "void config_set_color(Config*, const char*, Color*)"
  extern "const char* config_get_cstr(const Config*, const char*, " +
                                     "const char*, const char*)"
  extern "void config_set_cstr(const Config*, const char*, " +
                              "const char*, const char*)"
  extern "double config_get_num(const Config*, const char*, " +
                               "const char*, double)"
  extern "void config_set_num(const Config*, const char*, " +
                             "const char*, double)"
  extern "StrArray* config_get_cstr_list(const Config*, const char*, " +
                                        "const char*)"
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

    def get_color(key)
      color = GT::Color.malloc
      GT.config_get_colorptr(@config, color, key)
      color
    end

    def set_color(key, color)
      GT.config_set_color(@config, key, color)
    end

    def get_cstr(section, key)
      GT.config_get_cstr(@config, section, key, "undefined")
    end

    def set_cstr(section, key, value)
      GT.config_set_cstr(@config, section, key, value)
    end

    def get_num(section, key)
      GT.config_get_num(@config, section, key, -9999.9)
    end

    def set_num(section, key, number)
      GT.config_set_num(@config, section, key, number)
    end

    def get_cstr_list(section, key)
      strarray = GT.config_get_cstr_list(@config, section, key)
      if strarray then GT::StrArray.new(strarray) else nil end
    end

    def set_cstr_list(section, key, list)
      strarray = GT::StrArray.new()
      strarray.add_list(list)
      GT.config_set_cstr_list(@config, section, key, strarray.strarray)
    end
  end
end
