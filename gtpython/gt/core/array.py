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

class Array:
  def __init__(self, arr = None):
    from ctypes import sizeof, c_void_p
    if not arr:
      self.array = gtlib.gt_array_new(sizeof(c_void_p))
    else:
      self.array = arr
    self._as_parameter_ = self.array

  def __del__(self):
    gtlib.gt_array_delete(self.array)

  def from_param(cls, obj):
    if not isinstance(obj, Array):
      raise TypeError, "argument must be an Array"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def get(self, i):
    return gtlib.gt_array_get(self.array, i).contents

  def size(self):
    return gtlib.gt_array_size(self.array)

  def register(cls, gtlib):
    from ctypes import c_void_p, c_ulong, POINTER
    gtlib.gt_str_new_cstr.argtypes = [c_void_p, c_ulong]
    gtlib.gt_array_get.restype = POINTER(c_void_p)
    gtlib.gt_array_size.restype = c_ulong
  register = classmethod(register)
