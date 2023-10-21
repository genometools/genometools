#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Sascha Steinbiss <sascha@steinbiss.name>
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
from gt.core.error import Error, gterror


class TypeChecker(object):
    def __init__(self, *args):
        raise NotImplementedError(
            "Please call the constructor of a " + "TypeChecker implementation."
        )

    def __del__(self):
        try:
            gtlib.gt_type_checker_delete(self.tc)
        except AttributeError:
            pass

    def description(self):
        return str(gtlib.gt_type_checker_description(self.tc))

    def is_valid(self, t):
        return gtlib.gt_type_checker_is_valid(self.tc, gtbytes(t))

    def register(cls, gtlib):
        from ctypes import c_void_p, c_bool, c_char_p

        gtlib.gt_type_checker_is_valid.argtypes = [c_void_p, c_char_p]
        gtlib.gt_type_checker_is_valid.restype = c_bool
        gtlib.gt_type_checker_description.argtypes = [c_void_p]
        gtlib.gt_type_checker_description.restype = c_char_p
        gtlib.gt_type_checker_delete.argtypes = [c_void_p]

    register = classmethod(register)


class TypeCheckerOBO(TypeChecker):
    def __init__(self, filename):
        try:
            p = open(filename)
            tmp = p.readline()
            p.close()
        except:
            gterror("File " + filename + " not readable!")
        err = Error()
        self.tc = gtlib.gt_type_checker_obo_new(
            gtbytes(filename), err._as_parameter_)
        if self.tc == None:
            gterror(err)
        self._as_parameter_ = self.tc

    def from_param(cls, obj):
        if not isinstance(obj, TypeCheckerOBO):
            raise TypeError("argument must be a TypeCheckerOBO")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p

        gtlib.gt_type_checker_obo_new.argtypes = [c_char_p, c_void_p]
        gtlib.gt_type_checker_obo_new.restype = c_void_p

    register = classmethod(register)


class TypeCheckerBuiltin(TypeChecker):
    def __init__(self):
        self.tc = gtlib.gt_type_checker_builtin_new()
        self._as_parameter_ = self.tc

    def from_param(cls, obj):
        if not isinstance(obj, TypeCheckerOBO):
            raise TypeError("argument must be a TypeCheckerBuiltin")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_void_p

        gtlib.gt_type_checker_builtin_new.restype = c_void_p

    register = classmethod(register)
