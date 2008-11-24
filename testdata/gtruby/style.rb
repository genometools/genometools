#
# Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

# testing the Ruby bindings for the style object

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} style_file"
  STDERR.puts "Load style_file and test style bindings."
  exit(1)
end

stylefile = ARGV[0]

# create new style object
style = GT::Style.new()

# load style file
style.load_file(stylefile)

# clone style file
clone = style.clone
raise if not clone

# get color
color = style.get_color("exon", "fill")
raise if not color

# set color
color = GT::Color.malloc
color.red   = 0.3
color.green = 0.4
color.blue  = 0.3
style.set_color("exon", "fill", color)
color2 = style.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue

# unset color
style.unset("exon", "fill")
color2 = style.get_color("exon", "fill")
raise if not color2.nil?

# get undefined color
color = style.get_color("undefined", "undefined")
raise if not color.nil?

# get string
str = style.get_cstr("exon", "style")
raise if str != "box"

# set string
style.set_cstr("exon", "style", "line")
str = style.get_cstr("exon", "style")
raise if str != "line"

# unset string
style.unset("exon", "style")
str = style.get_cstr("exon", "style")
raise if not str.nil?

# get undefined string
str = style.get_cstr("undefined", "undefined")
raise if not str.nil?

# get number
num = style.get_num("format", "margins")
raise if num != 30

# set number
style.set_num("format", "margins", 20)
num = style.get_num("format", "margins")
raise if num != 20

# unset number
style.unset("format", "margins");
num = style.get_num("format", "margins")
raise if not num.nil?

#get undefined number
num = style.get_num("undefined", "undefined")
raise if not num.nil?

# get boolean
bool = style.get_bool("format", "show_grid")
raise if not bool

# get undefined boolean
bool = style.get_bool("undefined", "undefined")
raise if not bool.nil?

# set boolean
style.set_bool("format", "show_grid", false)
bool = style.get_bool("format", "show_grid")
raise if bool

# unset boolean
style.unset("format", "show_grid")
bool = style.get_bool("format", "show_grid")
raise if not bool.nil?

# serialise style to Lua code
style.set_num("format", "margins", 20)
style.set_bool("format", "show_grid", true)
color = GT::Color.malloc
color.red   = 0.3
color.green = 0.4
color.blue  = 0.3
style.set_color("exon", "fill", color)
luacode = style.to_str
raise if luacode.nil? or luacode.length == 0

# load style from Lua code
style.load_str(luacode)
num = style.get_num("format", "margins")
raise if num != 20
bool = style.get_bool("format", "show_grid")
raise if not bool
color2 = style.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue

# clone style from existing copy
style2 = style.clone
num = style2.get_num("format", "margins")
raise if num != 20
bool = style2.get_bool("format", "show_grid")
raise if not bool
color2 = style2.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue
style2.set_num("format", "margins", 30)
raise if style2.get_num("format", "margins") != 30
raise if style.get_num("format", "margins")\
           == style2.get_num("format", "margins")
