#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import gt
import os

op = os.path
datadir = op.abspath(op.join(op.dirname(__file__), "..", "..",
                     "testdata"))


class StreamTest(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "U89959_sas.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_pull(self):
        add_introns_stream = gt.AddIntronsStream(self.ins)
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(add_introns_stream, fi).pull()

        self.assert_('1877523' in fi.get_seqids())


if __name__ == "__main__":
    unittest.main()
