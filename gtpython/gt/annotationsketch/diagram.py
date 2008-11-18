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

from ctypes import CFUNCTYPE, c_char_p, c_void_p
from gt.dlload import gtlib
from gt.annotationsketch.block import Block
from gt.annotationsketch.canvas import Canvas
from gt.annotationsketch.feature_index import FeatureIndex
from gt.annotationsketch.style import Style
from gt.core.error import Error, gterror
from gt.core.range import Range

TrackSelectorFunc = CFUNCTYPE(c_char_p, c_void_p, c_void_p)

class Diagram:
  def __init__(self, feature_index, seqid, range, style):
    from ctypes import byref
    err = Error()
    if range.start > range.end:
      gterror("range.start > range.end")
    self.diagram = gtlib.gt_diagram_new(feature_index, seqid, byref(range), \
                                        style, err)
    if self.diagram == 0:
      gterror(err)
    self._as_parameter_ = self.diagram

  def __del__(self):
    try:
      gtlib.gt_diagram_delete(self.diagram)
    except AttributeError:
      pass

  def set_track_selector_func(self, func):
    def trackselector(block_ptr, data_ptr):
      b = Block(block_ptr)
      ret = func(b)
      if not ret:
        gterror("Track selector callback function must return a string!")
      return ret
    self.tsf_cb = TrackSelectorFunc(trackselector)
    self.tsf = trackselector
    gtlib.gt_diagram_set_track_selector_func(self.diagram, self.tsf_cb)

  def from_param(cls, obj):
    if not isinstance(obj, Diagram):
      raise TypeError, "argument must be a Diagram"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def register(cls, gtlib):
    from ctypes import c_char_p, c_void_p, POINTER
    gtlib.gt_diagram_new.restype = c_void_p
    gtlib.gt_diagram_new.argtypes = [c_void_p, c_char_p, POINTER(Range), \
                                     Style, Error]
    gtlib.gt_diagram_set_track_selector_func.argtypes = [c_void_p, \
                                                         TrackSelectorFunc]
  register = classmethod(register)
