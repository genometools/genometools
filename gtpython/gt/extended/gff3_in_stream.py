#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
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

from gt.bytes import gtbytes
from gt.dlload import gtlib
from gt.core.error import gterror
from gt.core.str_array import StrArray
from gt.extended.genome_stream import GenomeStream
from gt.extended.type_checker import TypeChecker


class GFF3InStream(GenomeStream):
    def __init__(self, filename):
        try:
            p = open(filename)
            tmp = p.readline()
            p.close()
        except:
            gterror("File " + filename + " not readable!")
        self.gs = gtlib.gt_gff3_in_stream_new_sorted(gtbytes(filename))
        self._as_parameter_ = self.gs

    def from_param(cls, obj):
        if not isinstance(obj, GFF3InStream):
            raise TypeError("argument must be a GFF3InStream")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def get_used_types(self):
        str_array_ptr = gtlib.gt_gff3_in_stream_get_used_types(self.gs)
        used_types = StrArray(str_array_ptr)
        return used_types.to_list()

    def check_id_attributes(self):
        gtlib.gt_gff3_in_stream_check_id_attributes(self.gs)

    def enable_tidy_mode(self):
        gtlib.gt_gff3_in_stream_enable_tidy_mode(self.gs)

    def enable_strict_mode(self):
        gtlib.gt_gff3_in_stream_enable_strict_mode(self.gs)

    def set_type_checker(self, tc):
        if not isinstance(tc, TypeChecker):
            raise TypeError("argument must be a TypeChecker")
        gtlib.gt_gff3_in_stream_set_type_checker(self.gs, tc._as_parameter_)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p

        gtlib.gt_gff3_in_stream_get_used_types.argtypes = [c_void_p]
        gtlib.gt_gff3_in_stream_new_sorted.argtypes = [c_char_p]
        gtlib.gt_gff3_in_stream_get_used_types.restype = c_void_p
        gtlib.gt_gff3_in_stream_new_sorted.restype = c_void_p
        gtlib.gt_gff3_in_stream_check_id_attributes.argtypes = [c_void_p]
        gtlib.gt_gff3_in_stream_enable_tidy_mode.argtypes = [c_void_p]
        gtlib.gt_gff3_in_stream_enable_strict_mode.argtypes = [c_void_p]
        gtlib.gt_gff3_in_stream_set_type_checker.argtypes = [
            c_void_p, c_void_p]

    register = classmethod(register)
