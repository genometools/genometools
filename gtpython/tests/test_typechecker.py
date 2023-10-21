#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import gt
import os

op = os.path
datadir = op.abspath(op.join(op.dirname(__file__), "..", "..", "testdata"))
gtdatadir = op.abspath(op.join(op.dirname(__file__), "..", "..", "gtdata"))


class TypeCheckerTest(unittest.TestCase):
    def test_builtin(self):
        tc = gt.extended.TypeCheckerBuiltin()
        self.assertTrue(tc.is_valid("gene"))
        self.assertFalse(tc.is_valid("nothing"))

    def test_obo_success(self):
        tc = gt.extended.TypeCheckerOBO(
            os.path.join(gtdatadir, "obo_files/sofa.obo"))
        self.assertTrue(tc.is_valid("gene"))
        self.assertFalse(tc.is_valid("nothing"))
        self.assertFalse(tc.is_valid(""))
        self.assertFalse(tc.is_valid(None))

    def test_obo_fail(self):
        self.assertRaises(
            gt.core.error.GTError, gt.extended.TypeCheckerOBO, "nonexistant"
        )
        self.assertRaises(
            gt.core.error.GTError,
            gt.extended.TypeCheckerOBO,
            os.path.join(datadir, "obo_files/corrupt_instance_stanza.obo"),
        )


if __name__ == "__main__":
    unittest.main()
