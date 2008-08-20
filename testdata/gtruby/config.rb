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

# testing the Ruby bindings for the Config object

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} config_file"
  STDERR.puts "Load config_file and test config bindings."
  exit(1)
end

configfile = ARGV[0]

# create new config object
config = GT::Config.new()

# load config file
config.load_file(configfile)

# get color
color = config.get_color("exon", "fill")
raise if not color

# set color
color = GT::Color.malloc
color.red   = 0.3
color.green = 0.4
color.blue  = 0.3
config.set_color("exon", "fill", color)
color2 = config.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue

# unset color
config.unset("exon", "fill")
color2 = config.get_color("exon", "fill")
raise if not color2.nil?

# get undefined color
color = config.get_color("undefined", "undefined")
raise if not color.nil?

# get string
str = config.get_cstr("exon", "style")
raise if str != "box"

# set string
config.set_cstr("exon", "style", "line")
str = config.get_cstr("exon", "style")
raise if str != "line"

# unset string
config.unset("exon", "style")
str = config.get_cstr("exon", "style")
raise if not str.nil?

# get undefined string
str = config.get_cstr("undefined", "undefined")
raise if not str.nil?

# get number
num = config.get_num("format", "margins")
raise if num != 30

# set number
config.set_num("format", "margins", 20)
num = config.get_num("format", "margins")
raise if num != 20

# unset number
config.unset("format", "margins");
num = config.get_num("format", "margins")
raise if not num.nil?

#get undefined number
num = config.get_num("undefined", "undefined")
raise if not num.nil?

# get boolean
bool = config.get_bool("format", "show_grid")
raise if not bool

# get undefined boolean
bool = config.get_bool("undefined", "undefined")
raise if not bool.nil?

# set boolean
config.set_bool("format", "show_grid", false)
bool = config.get_bool("format", "show_grid")
raise if bool

# unset boolean
config.unset("format", "show_grid")
bool = config.get_bool("format", "show_grid")
raise if not bool.nil?

# serialise Config to Lua code
config.set_num("format", "margins", 20)
config.set_bool("format", "show_grid", true)
color = GT::Color.malloc
color.red   = 0.3
color.green = 0.4
color.blue  = 0.3
config.set_color("exon", "fill", color)
luacode = config.to_str
raise if luacode.nil? or luacode.length == 0

# load config from Lua code
config.load_str(luacode)
num = config.get_num("format", "margins")
raise if num != 20
bool = config.get_bool("format", "show_grid")
raise if not bool
color2 = config.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue

# clone Config from existing copy
config2 = config.clone
num = config2.get_num("format", "margins")
raise if num != 20
bool = config2.get_bool("format", "show_grid")
raise if not bool
color2 = config2.get_color("exon", "fill")
raise if color2.red != color.red \
  and color2.green != color.green \
  and color2.blue != color.blue
config2.set_num("format", "margins", 30)
raise if config2.get_num("format", "margins") != 30
raise if config.get_num("format", "margins")\
           == config2.get_num("format", "margins")
