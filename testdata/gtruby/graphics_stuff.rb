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

require 'gtruby'

if ARGV.size != 2 then
  STDERR.puts "Usage: #{$0} curve_coords PNG_file"
  STDERR.puts "Tests drawing functions in the Graphics object."
  exit(1)
end

g = GT::GraphicsCairoSVG.new(700, 800)

black = [0,0,0,1].pack('dddd').to_ptr
black2 = [0,0,0,0.7].pack('dddd').to_ptr
red = [1,0,0,0.2].pack('dddd').to_ptr
red2 = [1,0,0,0.3].pack('dddd').to_ptr
grey1 = [0,0,0,0.3].pack('dddd').to_ptr
green = [0,1,0,0.2].pack('dddd').to_ptr
lastcol = [0,1,1,0.2].pack('dddd').to_ptr
blue = [0,0,1,0.6].pack('dddd').to_ptr

g.draw_text(300, 20, "draw_text")
g.draw_text_right(300, 40, "draw_text_left")
g.draw_text_centered(300, 60, "draw_text_centered")
g.draw_colored_text(300, 80, red, "draw_colored_text")
g.set_margins(20, 30)
raise if g.get_image_height() != 800
raise if g.get_image_width() != 700
raise if g.get_xmargins() != 20
raise if g.get_ymargins() != 30

g.draw_horizontal_line(150, 100, grey1, 300, 5.5)
g.draw_vertical_line(200, 70, red2, 100, 3.5)
g.draw_box(150, 120, 300, 30, green, GT::ARROW_LEFT,  20, 2, \
           black, false)
g.draw_box(150, 170, 300, 20, lastcol, GT::ARROW_RIGHT, 20, 1, \
           black, false)
g.draw_box(150, 210, 300, 20, red, GT::ARROW_BOTH,  20, 1, \
           black, false)
g.draw_box(300, 250,  20, 20, red, GT::ARROW_BOTH,  20, 1, \
           black, false)
g.draw_box(150, 290, 300, 20, red, GT::ARROW_NONE,  20, 1, \
           black, false)
g.draw_dashes(159, 340, 90, 20, GT::ARROW_RIGHT, 20, 1, black)
g.draw_dashes(359, 340, 60, 20, GT::ARROW_LEFT, 20, 1, black)
g.draw_dashes(59, 340, 60, 20, GT::ARROW_NONE, 20, 1, black)
g.draw_dashes(559, 340, 60, 20, GT::ARROW_BOTH, 20, 1, black)
g.draw_caret(159, 390, 90, 20, GT::ARROW_RIGHT, 20, 1, black)
g.draw_caret(359, 390, 60, 20, GT::ARROW_LEFT, 20, 1, black)
g.draw_caret(59, 390, 60, 20, GT::ARROW_NONE, 20, 1, black)
g.draw_caret(559, 390, 60, 20, GT::ARROW_BOTH, 20, 1, black)

data = []
File.open(ARGV[0]) do |file|
  file.each_line do |l|
    data.push(l.to_f)
  end
end
g.draw_curve_data(20, 430, blue, data, data.length, 0, 1, 40)

g.to_file(ARGV[1])
