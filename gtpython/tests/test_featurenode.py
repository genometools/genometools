
import unittest
from gt import FeatureNode, FeatureNodeIteratorDepthFirst

class FeatureNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode("test", "type", 100, 500, "+")

    def test_repr(self):
        self.assertEqual(str(self.fn),
            'FeatureNode(start=100, end=500, seqid="test")')


    def test_score(self):
        fn = self.fn
        self.assert_(not fn.score_is_defined())

        fn.set_score(2)
        self.assert_(fn.score_is_defined())
        self.assertEqual(2, fn.get_score())

        fn.unset_score()
        self.assert_(not fn.score_is_defined())

    def test_type(self):
        fn = self.fn
        self.assert_(not fn.has_type("foo"))
        self.assert_(fn.has_type("type"))

    def test_strand(self):
        fn = self.fn
        self.assertEqual(fn.get_strand(), "+")

    def test_seqid(self):
        fn = self.fn
        self.assertEqual(fn.seqid, "test")

    def test_start_end(self):
        fn = self.fn
        self.assertEqual(fn.start, 100)
        self.assertEqual(fn.end, 500)
       
    def test_attributes(self):
        fn = self.fn
        fn.add_attribute("test","testval")
        fn.add_attribute("test2","testval2")

        self.assert_("test" in fn.attribs)
        self.assert_("test2" in fn.attribs)

        nattrs = 0
        for tag, val in fn.each_attribute():
            self.assertEqual(val, fn.get_attribute(tag))
            nattrs += 1
        
        self.assertEqual(nattrs, 2)

class TestFeatureNodeChildren(unittest.TestCase):
    def setUp(self):
        self.fn  = FeatureNode("test", "type", 100, 500, "+")
        self.fn2 = FeatureNode("test", "type2", 200,300,"+")
        self.fn.add_child(self.fn2)

    def test_phase(self):
        fn = self.fn
        self.assertEqual(fn.get_phase(), 3)

        fn.set_phase(0)
        self.assertEqual(fn.get_phase(), 0)

    def test_fni(self):
        fn = self.fn
        fni = FeatureNodeIteratorDepthFirst(fn)
        num_features = 0
        tfn = fni.next()
        while tfn:
            tfn = fni.next()
            num_features += 1
        self.assertEqual(num_features, 2)



        fn3 = FeatureNode("test", "type3", 250,300,"+")
        fn.add_child(fn3)
        fni = FeatureNodeIteratorDepthFirst(fn)

        num_features = 0
        tfn = fni.next()
        while tfn:
            num_features += 1
            tfn = fni.next()
        self.assertEqual(num_features, 3)


class TestFeatureNodeProperties(unittest.TestCase):
    def setUp(self):
        self.fn  = FeatureNode("test", "type", 100, 500, "+")

    def test_strand(self):
        fn = self.fn
        self.assertEqual("+", fn.strand)
        fn.strand = "-"
        self.assertEqual("-", fn.strand)

    def test_score(self):
        fn = self.fn
        self.assert_(not fn.score_is_defined())
        fn.score = 2
        self.assert_(fn.score_is_defined())
        self.assertEqual(2, fn.get_score())
        self.assertEqual(2, fn.score)

        fn.set_score(4)

        self.assertEqual(2, fn.score)
    def test_range(self):
        fn = self.fn
        self.assertEqual((100, 500), fn.range)

if __name__ == "__main__":
    unittest.main()


