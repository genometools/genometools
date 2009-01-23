
import unittest
from gt import FeatureNode, FeatureNodeIteratorDepthFirst

class FeatureNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode.create("test", "type", 100, 500, "+")


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
        self.fn  = FeatureNode.create("test", "type", 100, 500, "+")
        self.fn2 = FeatureNode.create("test", "type2", 200,300,"+")
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



        fn3 = FeatureNode.create("test", "type3", 250,300,"+")
        fn.add_child(fn3)
        fni = FeatureNodeIteratorDepthFirst(fn)

        num_features = 0
        tfn = fni.next()
        while tfn:
            num_features += 1
            tfn = fni.next()
        self.assertEqual(num_features, 3)



if __name__ == "__main__":
    unittest.main()


