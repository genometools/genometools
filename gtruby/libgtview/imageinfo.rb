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
require 'dl/struct'
require 'gthelper'
require 'libgtext/genome_feature'
require 'libgtview/recmap'

module GT
  extend DL::Importable
  gtdlload "libgt"
  extern "ImageInfo* image_info_new()"
  extern "unsigned int image_info_get_height(ImageInfo*)"
  extern "unsigned long image_info_num_of_recmaps(ImageInfo*)"
  extern "void image_info_get_recmap_ptr(ImageInfo*, RecMap*, unsigned long)"
  extern "void image_info_delete(ImageInfo*)"

  class ImageInfo
    attr_reader :imageinfo
    def initialize()
      @imageinfo = GT.image_info_new()
      @imageinfo.free = GT::symbol("image_info_delete", "0P")
    end

    def get_height()
      GT.image_info_get_height(@imageinfo)
    end

    def num_of_recmaps()
      GT.image_info_num_of_recmaps(@imageinfo)
    end

    def each_hotspot()
      if @hotspots.nil? then
        @hotspots = []
        0.upto(self.num_of_recmaps()-1) do |i|
          rm = GT::RecMap.malloc
          GT.image_info_get_recmap_ptr(@imageinfo, rm, i)
          gf = GT::GenomeFeature.new(rm.gn, true)  #refcount only this GF!
          @hotspots.push([rm.nw_x.to_i, rm.nw_y.to_i, rm.se_x.to_i, rm.se_y.to_i, gf])
        end
        @hotspots.sort!{|hs1,hs2| hs1[2]-hs1[0]+1 <=> hs2[2]-hs2[0]+1}
      end
      @hotspots.each do |hs|
        yield hs[0],hs[1],hs[2],hs[3],hs[4]
      end
    end

    def to_ptr
      @imageinfo
    end
  end
end
