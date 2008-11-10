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
from gt.annotationsketch.rec_map import RecMap

class ImageInfo:
  def __init__(self):
    self.ii = gtlib.gt_image_info_new()
    self._as_parameter_ = self.ii
    self.hotspots = None

  def __del__(self):
    return gtlib.gt_image_info_delete(self.ii)

  def get_height(self):
    return gtlib.gt_image_info_get_height(self.ii)

  def num_of_rec_maps(self):
    return gtlib.gt_image_info_num_of_rec_maps(self.ii)

  def compare_hotspots(cls, hs1, hs2):
    if hs1[2]-hs1[0]+1 > hs2[2]-hs2[0]+1:
      return 1
    elif hs1[2] - hs1[0] + 1 == hs2[2] - hs2[0] + 1:
      return 0
    else:
      return -1
  compare_hotspots = classmethod(compare_hotspots)

  def each_hotspot(self):
    if not self.hotspots:
      self.hotspots = []
      for i in range(self.num_of_rec_maps()):
        rm = RecMap(gtlib.gt_image_info_get_rec_map(self.ii, i))
        self.hotspots.append([rm.get_northwest_x(), \
                              rm.get_northwest_y(), \
                              rm.get_southeast_x(), \
                              rm.get_southeast_y(), \
                              rm.get_genome_feature()])
      self.hotspots.sort(ImageInfo.compare_hotspots)
    for hs in self.hotspots:
      yield hs[0],hs[1],hs[2],hs[3],hs[4]

  def register(cls, gtlib):
    from ctypes import c_char_p
  register = classmethod(register)
