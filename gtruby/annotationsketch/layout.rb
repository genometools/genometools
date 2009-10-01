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

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtLayout*     gt_layout_new(GtDiagram*, unsigned int, GtStyle*, " + \
                                     "GtError*)"
  extern "unsigned long gt_layout_get_height(GtLayout*)"
  extern "int           gt_layout_sketch(GtLayout*, GtCanvas*, GtError*)"
  extern "void          gt_layout_set_track_ordering_func(GtLayout*, " + \
                                                          "void*, void*)"
  extern "void          gt_layout_delete(GtCanvas*)"

  class Layout
    def initialize(diagram, width, style)
      err = GT::Error.new()
      @layout = GT.gt_layout_new(diagram, width, style, err)
      if @layout.nil? then
        GT::gterror(err)
      end
      @layout.free = GT::symbol("gt_layout_delete", "0P")
    end

    def get_height()
      return GT.gt_layout_get_height(@layout)
    end

    def sketch(canvas)
      err = GT::Error.new()
      had_err = GT.gt_layout_sketch(@layout, canvas, err)
      if had_err < 0 then
        GT::gterror(err)
      end
    end

    def set_track_ordering_func(proc)
      @tof = DL.callback('ISSP') do |s1, s2, data_ptr|
               r = proc.call(s1, s2)
               if (!r.is_a?(Numeric)) then
                 GT::gterror("Track ordering callback must return a number!")
               end
               r.to_i
             end
      GT.gt_layout_set_track_ordering_func(@layout, @tof, GT::NULL)
    end

    def release_track_ordering_func
      DL.remove_callback(@tof)
    end

    def to_ptr
      @layout
    end
  end
end
