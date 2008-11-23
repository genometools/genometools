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
  extern "GtCustomTrack* gt_custom_track_script_wrapper_new(void*, void*, void*,
                                                            void*)"

  class CustomTrack
    def initialize()
      @get_height = DL.callback("L") do
                      puts "calling get_height"
                      self.get_height()
                    end
      @get_title = DL.callback("P") do
                      puts "calling get_title"
                      self.get_title().to_ptr
                    end
      @render    = DL.callback("IPLPPP") do |g, ypos, rng, sty, err|
                      puts "calling render"
                      self.render(GT::Graphics.new(g), ypos, rng,              \
                                  GT::Style.new(sty), GT::Error.new(err))
                    end
      @free    = DL.callback("0") do |g, ypos, rng, sty, err|
                      self.free()
                    end
      @ct = GT.gt_custom_track_script_wrapper_new(@render,                     \
                                                  @get_height,                 \
                                                  @get_title,                  \
                                                  @free)
      @ct.free = GT::symbol("gt_custom_track_delete", "0P")
    end

    def to_ptr
      @ct
    end
  end
end
