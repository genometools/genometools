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
from gt.annotationsketch.canvas import Canvas
from gt.annotationsketch.feature_index import FeatureIndex
from gt.annotationsketch.style import Style
from gt.core.error import Error, gterror
from gt.core.range import Range

class Diagram:
  def __init__(self, feature_index, seqid, range, style):
    from ctypes import byref
    if range.start > range.end:
      gterror("range.start > range.end")
    self.diagram = gtlib.gt_diagram_new(feature_index, seqid, byref(range), \
                                        style)
    self._as_parameter_ = self.diagram

  def __del__(self):
    try:
      gtlib.gt_diagram_delete(self.diagram)
    except AttributeError:
      pass

  def from_param(cls, obj):
    if not isinstance(obj, Diagram):
      raise TypeError, "argument must be a Diagram"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def sketch(self, canvas):
    return gtlib.gt_diagram_sketch(self.diagram, canvas)

  def register(cls, gtlib):
    from ctypes import c_char_p, c_void_p, POINTER
    gtlib.gt_diagram_new.restype = c_void_p
    gtlib.gt_diagram_new.argtypes = [FeatureIndex, c_char_p, POINTER(Range), \
                                     Style]
    gtlib.gt_diagram_sketch.argtypes = [c_void_p, Canvas]
  register = classmethod(register)
