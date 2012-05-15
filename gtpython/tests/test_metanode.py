#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import MetaNode, Str


class MetaNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = MetaNode.create_new("test", "data")
        self.fn2 = MetaNode.create_new(333, 444)

    def test_get_directive(self):
        self.assertEqual(self.fn.get_directive(), 'test')
        self.assertEqual(self.fn2.get_directive(), '333')

    def test_get_data(self):
        self.assertEqual(self.fn.get_data(), 'data')
        self.assertEqual(self.fn2.get_data(), '444')

if __name__ == "__main__":
    unittest.main()

