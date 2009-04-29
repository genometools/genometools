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


class TestFailedError(Exception):

    pass


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: " + (sys.argv)[0] +
                         " curve_coords PNG_file\n")
        sys.stderr.write("Tests drawing functions in the Graphics object.\n")
        sys.exit(1)
    g = GraphicsCairoSVG(700, 800)

    g.draw_text(300, 20, "draw_text")
    g.draw_text_right(300, 40, "draw_text_left")
    g.draw_text_centered(300, 60, "draw_text_centered")
    g.draw_colored_text(300, 80, Color(1, 0, 0, 0.7),
                        "draw_colored_text")
    g.set_margins(20, 30)
    if g.get_image_height() != 800:
        raise TestFailedError
    if g.get_image_width() != 700:
        raise TestFailedError
    if g.get_xmargins() != 20:
        raise TestFailedError
    if g.get_ymargins() != 30:
        raise TestFailedError

    g.draw_horizontal_line(150, 100, Color(0, 0, 0, 0.3), 300, 5.5)
    g.draw_vertical_line(200, 70, Color(1, 0, 0, 0.3), 100, 3.5)
    g.draw_box(
        150,
        120,
        300,
        30,
        Color(0, 1, 0, 0.2),
        ARROW_LEFT,
        20,
        2,
        Color(0, 0, 0, 1),
        False,
        )
    g.draw_box(
        150,
        170,
        300,
        20,
        Color(0, 1, 1, 0.2),
        ARROW_RIGHT,
        20,
        1,
        Color(0, 0, 0, 1),
        False,
        )
    g.draw_box(
        150,
        210,
        300,
        20,
        Color(1, 0, 0, 0.2),
        ARROW_BOTH,
        20,
        1,
        Color(0, 0, 0, 1),
        False,
        )
    g.draw_box(
        300,
        250,
        20,
        20,
        Color(1, 0, 0, 0.2),
        ARROW_BOTH,
        20,
        1,
        Color(0, 0, 0, 1),
        False,
        )
    g.draw_box(
        150,
        290,
        300,
        20,
        Color(1, 0, 0, 0.2),
        ARROW_NONE,
        20,
        1,
        Color(0, 0, 0, 1),
        False,
        )
    g.draw_dashes(159, 340, 90, 20, ARROW_RIGHT, 20, 1, Color(0, 0, 0, 1))
    g.draw_dashes(359, 340, 60, 20, ARROW_LEFT, 20, 1, Color(0, 0, 0, 1))
    g.draw_dashes(59, 340, 60, 20, ARROW_NONE, 20, 1, Color(0, 0, 0, 1))
    g.draw_dashes(559, 340, 60, 20, ARROW_BOTH, 20, 1, Color(0, 0, 0, 1))
    g.draw_caret(159, 390, 90, 20, ARROW_RIGHT, 20, 1, Color(0, 0, 0, 1))
    g.draw_caret(359, 390, 60, 20, ARROW_LEFT, 20, 1, Color(0, 0, 0, 1))
    g.draw_caret(59, 390, 60, 20, ARROW_NONE, 20, 1, Color(0, 0, 0, 1))
    g.draw_caret(559, 390, 60, 20, ARROW_BOTH, 20, 1, Color(0, 0, 0, 1))

    from random import random
    data = []
    ndata = 0
    f = open((sys.argv)[1])
    for line in f:
        data.append(float(line))
        ndata += 1
    g.draw_curve_data(20, 430, Color(0, 0, 1, .6), data, ndata, Range(0,
                      1), 40)

    g.to_file((sys.argv)[2])
