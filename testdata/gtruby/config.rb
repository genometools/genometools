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
color = config.get_color("exon")
raise if not color

# set color
color = GT::Color.malloc
color.red   = 0.3
color.green = 0.4
color.blue  = 0.3
config.set_color("exon", color)

# get string
str = config.get_cstr("feature_styles", "exon")
raise if str != "box"

# set string
config.set_cstr("feature_styles", "exon", "line")

# get number
num = config.get_num("format", "margins")
raise if num != 30

# set number
config.set_num("format", "margins", 20.0)

# get string list
strarray = config.get_cstr_list("collapse", "to_parent")
list = strarray.to_a # convert StrArray C object to Ruby array of strings
puts list

# set string list
list = [ "mRNA", "gene" ]
config.set_cstr_list("dominate", "exon", list)
