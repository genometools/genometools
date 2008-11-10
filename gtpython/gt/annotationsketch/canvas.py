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

from gt.dlload import gtlib
from gt.core.error import Error, gterror
from gt.core.str import Str

class Canvas:
  def __init__(self, *args):
    raise NotImplementedError, "Please call the constructor of a " \
                               "Canvas implementation."

  def __del__(self):
    gtlib.gt_canvas_delete(self.canvas)

class CanvasCairoFile(Canvas):
  def __init__(self, style, width, ii):
    self.canvas = gtlib.gt_canvas_cairo_file_new(style, 1, width, ii)
    self._as_parameter_ = self.canvas

  def to_file(self, filename):
    err = Error()
    rval = gtlib.gt_canvas_cairo_file_to_file(self.canvas, filename, err)
    if rval != 0:
      gterror(err)

  def to_stream(self):
    from ctypes import string_at
    str = Str(None)
    gtlib.gt_canvas_cairo_file_to_stream(self.canvas, str)
    return string_at(str.get_mem(), str.length())

  def register(cls, gtlib):
    from ctypes import c_char_p, c_void_p
    gtlib.gt_canvas_cairo_file_to_stream.restype  = c_char_p
    gtlib.gt_canvas_cairo_file_new.restype  = c_void_p
  register = classmethod(register)
