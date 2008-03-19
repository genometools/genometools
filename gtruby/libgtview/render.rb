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

require 'dl/import'
require 'gthelper'
require 'libgtcore/str'

module GT
  extend DL::Importable
  gtdlload "libgt"
  extern "Render* render_new(Config*)"
  extern "int render_to_png(Render*, Diagram*, const char*, unsigned int, " +
                           "Error*)"
  extern "void render_to_png_stream(Render*, Diagram*, Str*, unsigned int)"
  extern "void render_delete(Render*)"

  class Render
    attr_reader :render
    def initialize(config)
      @render = GT.render_new(config.config)
      @render.free = GT::symbol("render_delete", "0P")
    end

    def to_png(diagram, filename, width = 800)
      err = GT::Error.new()
      rval = GT.render_to_png(@render, diagram.diagram, filename, width,
                              err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def to_png_stream(diagram, width = 800)
      str = GT::Str.new(nil)
      GT.render_to_png_stream(@render, diagram.diagram, str.to_ptr, width)
      str.get_mem.to_s(str.length)
    end
  end
end
