#!/usr/bin/python
# -*- coding: utf-8 -*-
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

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
import sys


class TestFailedException(Exception):

    pass


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: " + (sys.argv)[0] + " style_file")
        sys.stderr.write("Load style_file and test style bindings.")
        sys.exit(1)

    stylefile = (sys.argv)[1]

  # create new style object

    style = Style()

  # load style file

    style.load_file(stylefile)

  # clone style file

    clone = style.clone()
    if not clone:
        raise TestFailedException

  # get color

    color = style.get_color("exon", "fill")
    if not color:
        raise TestFailedException

  # set color

    color = Color(0.3, 0.4, 0.3)
    style.set_color("exon", "fill", color)
    color2 = style.get_color("exon", "fill")
    if color2.red != color.red and color2.green != color.green and \
        color2.blue != color.blue:
        raise TestFailedException

  # unset color

    style.unset("exon", "fill")
    color2 = style.get_color("exon", "fill")
    if color2:
        raise TestFailedException

  # get undefined color

    color = style.get_color("undefined", "undefined")
    if color:
        raise TestFailedException

  # get string

    string = style.get_cstr("exon", "style")
    if string != "box":
        raise TestFailedException

  # set string

    style.set_cstr("exon", "style", "line")
    string = style.get_cstr("exon", "style")
    if string != "line":
        raise TestFailedException

  # unset string

    style.unset("exon", "style")
    string = style.get_cstr("exon", "style")
    if string:
        raise TestFailedException

  # get undefined string

    string = style.get_cstr("undefined", "undefined")
    if string:
        raise TestFailedException

  # get number

    num = style.get_num("format", "margins")
    if num != 30:
        raise TestFailedException

  # set number

    style.set_num("format", "margins", 20)
    num = style.get_num("format", "margins")
    if num != 20:
        raise TestFailedException

  # unset number

    style.unset("format", "margins")
    None
    num = style.get_num("format", "margins")
    if num:
        raise TestFailedException

  #get undefined number

    num = style.get_num("undefined", "undefined")
    if num:
        raise TestFailedException

  # get boolean

    bool = style.get_bool("format", "show_grid")
    if bool == False:
        raise TestFailedException

  # get undefined boolean

    bool = style.get_bool("undefined", "undefined")
    if not bool == None:
        raise TestFailedException

  # set boolean

    style.set_bool("format", "show_grid", False)
    bool = style.get_bool("format", "show_grid")
    if bool:
        raise TestFailedException

  # unset boolean

    style.unset("format", "show_grid")
    bool = style.get_bool("format", "show_grid")
    if not bool == None:
        raise TestFailedException

  # serialise style to Lua code

    style.set_num("format", "margins", 20)
    style.set_bool("format", "show_grid", True)
    color = Color(0.3, 0.4, 0.3)
    style.set_color("exon", "fill", color)
    luacode = style.to_str()
    if luacode == None or len(luacode) == 0:
        raise TestFailedException

  # load style from Lua code

    style.load_str(luacode)
    num = style.get_num("format", "margins")
    if num != 20:
        raise TestFailedException
    bool = style.get_bool("format", "show_grid")
    if not bool:
        raise TestFailedException
    color2 = style.get_color("exon", "fill")
    if color2.red != color.red and color2.green != color.green and \
        color2.blue != color.blue:
        raise TestFailedException

  # clone style from existing copy

    style2 = style.clone()
    num = style2.get_num("format", "margins")
    if num != 20:
        raise TestFailedException
    bool = style2.get_bool("format", "show_grid")
    if not bool:
        raise TestFailedException
    color2 = style2.get_color("exon", "fill")
    if color2.red != color.red and color2.green != color.green and \
        color2.blue != color.blue:
        raise TestFailedException
    style2.set_num("format", "margins", 30)
    if style2.get_num("format", "margins") != 30:
        raise TestFailedException
    if style.get_num("format", "margins") == style2.get_num("format",
            "margins"):
        raise TestFailedException
