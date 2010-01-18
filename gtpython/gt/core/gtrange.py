#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg
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
from gt.core.error import gterror
from ctypes import Structure, c_ulong


class Range(Structure):

    _fields_ = [("w_start", c_ulong), ("w_end", c_ulong)]

    def __init__(self, start = 0, end = 0):
        if start > end or start < 0 or end < 0:
            gterror("range error: start > end!")
        super(Range, self).__init__(start, end)

    @property
    def start(self):
        return self.w_start

    @start.setter
    def start(self, val):
        if val > self.end or not val > 0:
            gterror("Invalid range start component: %d" % val)
        self.w_start = val

    @property
    def end(self):
        return self.w_end

    @end.setter
    def end(self, val):
        if val < self.start or not val > 0:
            gterror("Invalid range end component: %d" % val)
        self.w_end = val
