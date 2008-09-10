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

require 'dl/import'
require 'gthelper'
require 'libgtcore/str'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GT_Canvas* gt_canvas_cairo_file_new(GT_Style*, int, unsigned int, "+
                                          "GT_ImageInfo*)"
  extern "int gt_canvas_cairo_file_to_file(GT_CanvasCairoFile*, const char*, " +
                                       "GT_Error*)"
  extern "int gt_canvas_cairo_file_to_stream(GT_CanvasCairoFile*, GT_Str*)"
  extern "void gt_canvas_delete(GT_Canvas*)"

  class Canvas
    def initialize(style, width, ii)
      @canvas = GT.gt_canvas_cairo_file_new(style.config, 1, width, ii.to_ptr)
      @canvas.free = GT::symbol("gt_canvas_delete", "0P")
    end

    def to_file(filename)
      err = GT::Error.new()
      rval = GT.gt_canvas_cairo_file_to_file(@canvas, filename, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def to_stream()
      str = GT::Str.new(nil)
      GT.gt_canvas_cairo_file_to_stream(@canvas, str.to_ptr)
      str.get_mem.to_s(str.length)
    end

    def to_ptr
      @canvas
    end
  end
end
