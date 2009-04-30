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
import re


class TestFailedError(Exception):

    pass


if __name__ == "__main__":
    fn = FeatureNode.create_new("test", "type", 100, 500, "+")
    fn.add_attribute("test", "testval")
    fn.add_attribute("test2", "testval2")
    if fn.score_is_defined():
        raise TestFailedError
    fn.set_score(2)
    if not fn.score_is_defined():
        raise TestFailedError
    if fn.get_score() != 2:
        raise TestFailedError
    fn.unset_score()
    if fn.score_is_defined():
        raise TestFailedError
    if fn.has_type("foo"):
        raise TestFailedError
    if not fn.has_type("type"):
        raise TestFailedError
    if not fn.get_type() == "type":
        raise TestFailedError
    fn.set_type("foo")
    if not fn.has_type("foo"):
        raise TestFailedError
    if not fn.get_type() == "foo":
        raise TestFailedError
    if not fn.get_strand() == "+":
        raise TestFailedError
    fn2 = FeatureNode.create_new("test", "type2", 200, 300, "+")
    fn.add_child(fn2)
    num_attrs = 0
    for (tag, val) in fn.each_attribute():
        if not val == fn.get_attribute(tag):
            raise TestFailedError
        num_attrs += 1
    if not num_attrs == 2:
        raise TestFailedError
    if fn.get_phase() != 3:  #undefined
        raise TestFailedError
    fn.set_phase(0)
    if fn.get_phase() != 0:  #zero
        raise TestFailedError

    fni = FeatureNodeIteratorDepthFirst(fn)
    num_features = 0
    tfn = fni.next()
    while tfn:
        num_features += 1
        tfn = fni.next()
    if num_features != 2:
        raise TestFailedError

    fn3 = FeatureNode.create_new("test", "type3", 250, 300, "+")
    fn.add_child(fn3)
    fni = FeatureNodeIteratorDepthFirst(fn)
    num_features = 0
    tfn = fni.next()
    while tfn:
        num_features += 1
        tfn = fni.next()
    if num_features != 3:
        raise TestFailedError
