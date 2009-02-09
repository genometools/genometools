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
from gt.core.gtstr import Str
from gt.core.str_array import StrArray
from gt.extended.genome_node import GenomeNode
from gt.props import cachedproperty

class FeatureNode(GenomeNode):

  @classmethod
  def create_from_ptr(cls, node_ptr, newref=False):
    fn = GenomeNode.__new__(cls, node_ptr, newref)
    fn.gn = node_ptr
    fn._as_parameter_ = fn.gn
    fn.attribs = {}
    fn.update_attrs()
    return fn

  def __init__(self, seqid, type, start, end, strand):
    from gt.extended.strand import strandchars
    if not strand in strandchars:
      gterror("Invalid strand '%s' -- must be one of %s" \
                 % (strand, strandchars))
    s = Str(seqid)
    newfn = gtlib.gt_feature_node_new(s, type, start, end, \
                                      strandchars.index(strand))
    super(FeatureNode, self).__init__(newfn, True)
    self.attribs = {}
    self.update_attrs()


  def update_attrs(self):
    def py_collect_func(tag, val, data):
      self.attribs[tag] = val
    self.collect_func = CollectFunc(py_collect_func)
    gtlib.gt_feature_node_foreach_attribute(self.gn, \
                                            self.collect_func, \
                                            None)

  def add_child(self, node):
    gtlib.gt_feature_node_add_child(self.gn, node)

  def from_param(cls, obj):
    if not isinstance(obj, FeatureNode):
      raise TypeError, "argument must be a FeatureNode"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def get_source(self):
    return gtlib.gt_feature_node_get_source(self.gn)

  def set_source(self, source):
    s = Str(source)
    gtlib.gt_feature_node_set_source(self.gn, s)
  source = cachedproperty(get_source, set_source)

  def get_type(self):
    return gtlib.gt_feature_node_get_type(self.gn)
  type = cachedproperty(get_type, None)

  def has_type(self, type):
    return (gtlib.gt_feature_node_has_type(self.gn, type) == 1)

  def set_strand(self, strand):
    from gt.extended.strand import strandchars
    if not strand in strandchars:
      gterror("Invalid strand '%s' -- must be one of %s" \
                 % (strand, strandchars))
    gtlib.gt_feature_node_set_strand(self.gn, strandchars.index(strand))

  def get_strand(self):
    from gt.extended.strand import strandchars
    return strandchars[gtlib.gt_feature_node_get_strand(self.gn)]
  strand = cachedproperty(get_strand, set_strand)

  def get_phase(self):
    return gtlib.gt_feature_node_get_phase(self.gn)

  def set_phase(self, phase):
    return gtlib.gt_feature_node_set_phase(self.gn, phase)
  phase = cachedproperty(get_phase, set_phase)

  def score_is_defined(self):
    return (gtlib.gt_feature_node_score_is_defined(self.gn) == 1)

  def get_score(self):
    if gtlib.gt_feature_node_score_is_defined(self.gn) == 1:
      return gtlib.gt_feature_node_get_score(self.gn)
    else:
      return None

  def set_score(self, score):
    gtlib.gt_feature_node_set_score(self.gn, score)

  def unset_score(self):
    gtlib.gt_feature_node_unset_score(self.gn)
  score = cachedproperty(get_score, set_score, unset_score)

  def get_attribute(self, attrib):
    return self.attribs[attrib]

  def add_attribute(self, attrib, value):
    if attrib == "" or value == "":
      gterror("attribute keys or values must not be empty!")
    gtlib.gt_feature_node_add_attribute(self.gn, attrib, value)
    self.update_attrs()

  def each_attribute(self):
    for tag, val in self.attribs.iteritems():
      yield (tag, val)

  def register(cls, gtlib):
    from ctypes import c_char_p, c_float, c_int, c_int, c_void_p, c_ulong, \
                       c_float
    gtlib.gt_feature_node_new.restype = c_void_p
    gtlib.gt_feature_node_new.argtypes = [Str, c_char_p, c_ulong, c_ulong,\
                                          c_int]
    gtlib.gt_feature_node_add_child.argtypes = [c_void_p, FeatureNode]
    gtlib.gt_feature_node_set_source.argtypes = [c_void_p, Str]
    gtlib.gt_feature_node_get_source.restype = c_char_p
    gtlib.gt_feature_node_get_source.argtypes = [c_void_p]
    gtlib.gt_feature_node_get_type.restype = c_char_p
    gtlib.gt_feature_node_get_type.argtypes = [c_void_p]
    gtlib.gt_feature_node_has_type.restype = c_int
    gtlib.gt_feature_node_has_type.argtypes = [c_void_p, c_char_p]
    gtlib.gt_feature_node_get_score.restype = c_float
    gtlib.gt_feature_node_get_score.argtypes = [c_void_p]
    gtlib.gt_feature_node_set_score.argtypes = [c_void_p, c_float]
    gtlib.gt_feature_node_get_phase.restype = c_int
    gtlib.gt_feature_node_get_phase.argtypes = [c_void_p]
    gtlib.gt_feature_node_set_phase.argtypes = [c_void_p, c_int]
    gtlib.gt_feature_node_score_is_defined.restype = c_int
    gtlib.gt_feature_node_score_is_defined.argtypes = [c_void_p]
    gtlib.gt_feature_node_get_strand.restype = c_int
    gtlib.gt_feature_node_get_strand.argtypes = [c_void_p]
    gtlib.gt_feature_node_set_strand.argtypes = [c_void_p, c_int]
    gtlib.gt_feature_node_unset_score.argtypes = [c_void_p]
    gtlib.gt_feature_node_add_attribute.argtypes = [c_void_p, c_char_p, \
                                                    c_char_p]
  register = classmethod(register)

class FeatureNodeIterator(object):
  def next(self):
    ret = gtlib.gt_feature_node_iterator_next(self.i)
    if ret != None:
      return FeatureNode.create_from_ptr(ret)
    return ret

  def __del__(self):
    try:
      gtlib.gt_feature_node_iterator_delete(self.i)
    except AttributeError:
      pass

  def register(cls, gtlib):
    from ctypes import c_void_p
    gtlib.gt_feature_node_iterator_new.restype = c_void_p
    gtlib.gt_feature_node_iterator_new.argtypes = [FeatureNode]
    gtlib.gt_feature_node_iterator_new_direct.restype = c_void_p
    gtlib.gt_feature_node_iterator_new_direct.argtypes = [FeatureNode]
    gtlib.gt_feature_node_iterator_next.restype = c_void_p
  register = classmethod(register)

class FeatureNodeIteratorDepthFirst(FeatureNodeIterator):
  def __init__(self, node):
    self.i = gtlib.gt_feature_node_iterator_new(node)
    self._as_parameter_ = self.i

class FeatureNodeIteratorDirect(FeatureNodeIterator):
  def __init__(self, node):
    self.i = gtlib.gt_feature_node_iterator_new_direct(node)
    self._as_parameter_ = self.i
