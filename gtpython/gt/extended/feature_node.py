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

from gt.dlload import gtlib, CollectFunc
from gt.core.error import Error, gterror
from gt.core.str_array import StrArray
from gt.extended.genome_node import GenomeNode

def py_collect_func(tag, val, data):
  gtlib.gt_str_array_add_cstr(data, tag)
  gtlib.gt_str_array_add_cstr(data, val)
collect_func = CollectFunc(py_collect_func)

class FeatureNode(GenomeNode):
  def __init__(self, node_ptr, single = False):
    super(FeatureNode, self).__init__(node_ptr, single)
    a = StrArray()
    gtlib.gt_feature_node_foreach_attribute(self.gn, \
                                            collect_func, \
                                            a)
    attribs_a = a.to_list()
    self.attribs = {}
    while len(attribs_a) > 0:
      tagname = attribs_a.pop(0)
      tagval = attribs_a.pop(0)
      self.attribs[tagname] = tagval

  def from_param(cls, obj):
    if not isinstance(obj, FeatureNode):
      raise TypeError, "argument must be a FeatureNode"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def get_type(self):
    return gtlib.gt_feature_node_get_type(self.gn)

  def get_strand(self):
    return gtlib.gt_feature_node_get_strand(self.gn)

  def get_phase(self):
    return gtlib.gt_feature_node_get_phase(self.gn)

  def get_score(self):
    if gtlib.gt_feature_node_score_is_defined(self.gn):
      return gtlib.gt_feature_node_get_score(self.gn)
    else:
      return None

  def get_attribute(self, attrib):
    return self.attribs[attrib]

  def each_attribute(self):
    for tag, val in self.attribs:
      yield tag, val

  def register(cls, gtlib):
    from ctypes import c_char_p, c_float, c_bool, c_int, c_void_p
    gtlib.gt_feature_node_get_type.restype = c_char_p
    gtlib.gt_feature_node_get_type.argtypes = [c_void_p]
    gtlib.gt_feature_node_get_score.restype = c_float
    gtlib.gt_feature_node_get_score.argtypes = [c_void_p]
    gtlib.gt_feature_node_get_phase.restype = c_int
    gtlib.gt_feature_node_get_phase.argtypes = [c_void_p]
    gtlib.gt_feature_node_get_strand.restype = c_int
    gtlib.gt_feature_node_get_strand.argtypes = [c_void_p]
    gtlib.gt_feature_node_score_is_defined.restype = c_bool
    gtlib.gt_feature_node_score_is_defined.argtypes = [c_void_p]
  register = classmethod(register)
