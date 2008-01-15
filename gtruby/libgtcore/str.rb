#
# Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

require 'dl/import'

module GT
  extend DL::Importable
  dlload "libgtcore.so"
  extern "Str* str_new()"
  extern "void* str_get_mem(const Str*)"
  extern "unsigned long str_length(const Str*)"
  extern "void str_delete(Str*)"

  class Str
    def initialize
      @str = GT.str_new()
      @str.free = GT::symbol("str_delete", "0P")
    end

    def get_mem
      GT.str_get_mem(@str)
    end

    def to_ptr
      @str
    end

    def length
      GT.str_length(@str)
    end
  end
end
