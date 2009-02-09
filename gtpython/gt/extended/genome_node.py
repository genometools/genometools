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
from gt.extended.gff3_visitor import GFF3Visitor
from gt.core.gtstr import Str
from gt.props import cachedproperty

class GenomeNode(object):
  def __init__(self, node_ptr, newref=False):
    if node_ptr == 0 or node_ptr == None:
      gterror("GenomeNode pointer cannot be NULL (was: " + str(node_ptr) + ")")
    if newref:
      self.gn = gtlib.gt_genome_node_ref(node_ptr)
    else:
      self.gn = node_ptr
    self._as_parameter_ = self.gn

  @classmethod
  def create_from_ptr(cls, node_ptr, newref=False):
    n = cls(node_ptr, newref)
    n.gn = newref and gtlib.genome_node_ref(node_ptr) or node_ptr
    n._as_parameter_ = n.gn
    return n

  def __repr__(self):
      c = self.__class__.__name__
      return "%s(start=%i, end=%i, seqid=\"%s\")" % \
              (c, self.start, self.end, self.seqid)

  def __del__(self):
    try:
      gtlib.gt_genome_node_delete(self.gn)
    except AttributeError:
      pass

  def from_param(cls, obj):
    if not isinstance(obj, GenomeNode):
      raise TypeError, "argument must be a GenomeNode"
    return obj._as_parameter_
  from_param = classmethod(from_param)

  def get_range(self):
    return (gtlib.gt_genome_node_get_start(self.gn), \
            gtlib.gt_genome_node_get_end(self.gn))
  range = property(get_range)

  def get_seqid(self):
    return str(Str(gtlib.gt_genome_node_get_seqid(self.gn)))
  seqid = cachedproperty(get_seqid)

  def get_start(self):
    return gtlib.gt_genome_node_get_start(self.gn)
  start = cachedproperty(get_start)

  def get_end(self):
    return gtlib.gt_genome_node_get_end(self.gn)
  end = cachedproperty(get_end)

  def get_filename(self):
    return gtlib.gt_genome_node_get_filename(self.gn)
  filename = property(get_filename)

  def get_line_number(self):
    return gtlib.gt_genome_node_get_line_number(self.gn)
  line_number = property(get_line_number)

  def accept(self, visitor):
    err = Error()
    rval = gtlib.gt_genome_node_accept(self.gn, visitor, err)
    if rval != 0:
      gterror(err)

  def register(cls, gtlib):
    from ctypes import c_char_p, c_ulong, c_int, c_void_p, c_uint
    gtlib.gt_genome_node_get_filename.restype = c_char_p
    gtlib.gt_genome_node_get_filename.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_start.restype = c_ulong
    gtlib.gt_genome_node_get_start.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_end.restype = c_ulong
    gtlib.gt_genome_node_get_end.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_seqid.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_seqid.restype = Str
    gtlib.gt_genome_node_get_filename.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_filename.restype = c_char_p
    gtlib.gt_genome_node_get_line_number.argtypes = [c_void_p]
    gtlib.gt_genome_node_get_line_number.restype = c_uint
    gtlib.gt_genome_node_accept.restype = c_int
    gtlib.gt_genome_node_accept.argtypes = [c_void_p, GFF3Visitor, Error]
  register = classmethod(register)
