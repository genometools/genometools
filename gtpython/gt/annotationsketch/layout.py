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

from gt.dlload import gtlib
from gt.annotationsketch.canvas import Canvas
from gt.annotationsketch.diagram import Diagram
from gt.annotationsketch.style import Style
from gt.core.error import Error, gterror


class Layout:

    def __init__(self, diagram, width, style):
        err = Error()
        self.layout = gtlib.gt_layout_new(diagram, width, style, err)
        if err.is_set():
            gterror(err)
        self._as_parameter_ = self.layout

    def __del__(self):
        try:
            gtlib.gt_layout_delete(self.layout)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, Layout):
            raise TypeError, "argument must be a Layout"
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def sketch(self, canvas):
        err = Error()
        had_err = gtlib.gt_layout_sketch(self.layout, canvas, err)
        if had_err < 0:
            gterror(err)

    def get_height(self):
        return gtlib.gt_layout_get_height(self.layout)

    def register(cls, gtlib):
        from ctypes import c_ulong, c_void_p, c_int
        gtlib.gt_layout_new.restype = c_void_p
        gtlib.gt_layout_new.argtypes = [Diagram, c_ulong, Style]
        gtlib.gt_layout_sketch.restype = c_int
        gtlib.gt_layout_sketch.argtypes = [c_void_p, Canvas, Error]
        gtlib.gt_layout_get_height.restype = c_ulong
        gtlib.gt_layout_get_height.argtypes = [c_void_p]

    register = classmethod(register)


