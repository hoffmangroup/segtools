import inspect
import sys
import unittest

from numpy import array
from tempfile import NamedTemporaryFile

from .common import load_features, load_segmentation
from .overlap import calc_overlap

def data2bed(lines):
    # Input lines should each be a tuple of (chrom, start, end, label)
    strings = []
    for line in lines:
        strings.append("chr%s\t%d\t%d\t%s" % line)
    return "\n".join(strings)

def data2gff(lines):
    # Input lines should each be a tuple of (chrom, start, end, label)
    strings = []
    for line in lines:
        chr, start, end, label = line
        strings.append("chr%s\t%s\t%s\t%d\t%d" % \
                           (chr, "source", label, start+1, end))
    return "\n".join(strings)

class OverlapTester(unittest.TestCase):
    def init(self):
        pass

    def setUp(self):
        self.kwargs = {}
        self.init()

        # Set self.features from self.subject, self.query
        self.features = []
        self._open_files = []
        for type, data in [self.subject, self.query]:
            new_file = NamedTemporaryFile()
            self._open_files.append(new_file)
            features = None
            if type == "bed":
                new_file.write(data2bed(data))
                new_file.flush()
                features = load_segmentation(new_file.name, verbose=False)
            elif type == "gff":
                new_file.write(data2gff(data))
                new_file.flush()
                features = load_features(new_file.name, verbose=False)
            elif type == "gtf":
                new_file.write(data2gff(data))
                new_file.flush()
                features = load_features(new_file.name, gtf=True, sort=True,
                                         verbose=False)

            if features:
                self.features.append(features)

    def tearDown(self):
        for file in self._open_files:
            file.close()

    def test(self):
        stats = calc_overlap(*self.features, **self.kwargs)
        self.assertValuesEqual(stats, self.answer)

    def assertArraysEqual(self, arr1, arr2):
        not_equal = (arr1 != arr2)
        if not_equal.any():
            self.fail("%r != %r" % (arr1, arr2))

    def assertValuesEqual(self, observed, expected):
        # counts, totals, nones
        for val1, val2 in zip(observed, expected):
            self.assertArraysEqual(val1, array(val2))

class TestPerfectOverlap(OverlapTester):
    def init(self):
        self.subject = ("bed", [(1, 0, 50, 1),
                                (1, 50, 100, 2)])
        self.query = ("bed", [(1, 0, 50, 1),
                              (1, 50, 100, 2)])
        self.answer = ([[1, 0], [0, 1]],
                       (1, 1),
                       (0, 0))

class TestNoOverlap(OverlapTester):
    def init(self):
        self.subject = ("bed", [(1, 10, 20, 1)])
        self.query = ("bed", [(1, 20, 35, 6)])
        self.answer = ([0],
                       (1),
                       (1))

class TestSimpleOverlap(OverlapTester):
    def init(self):
        self.subject = ("bed", [(1, 11, 20, 1),
                                (1, 25, 35, 2),
                                (1, 37, 40, 3)])
        self.query = ("gff", [(1, 15, 37, 2)])
        self.answer = ([[1], [1], [0]],
                       (1, 1, 1),
                       (0, 0, 1))

class TestStackedFeatures(OverlapTester):
    def init(self):
        self.subject = ("bed", [(1, 0, 10, 1)])
        self.query = ("gff", [(1, 2, 9, 2),
                              (1, 3, 12, 3),
                              (1, 4, 7, 4),
                              (1, 5, 6, 5)])
        self.answer = ([1, 1, 1, 1],
                       (1),
                       (0))

class TestOverlappingFeatures(OverlapTester):
    def init(self):
        self.subject = ("bed", [(1, 0, 10, 1),
                                (1, 10, 20, 2),
                                (1, 20, 30, 1)])
        self.query = ("gff", [(1, 0, 30, 1),
                              (1, 2, 15, 1),
                              (1, 5, 20, 2),
                              (1, 10, 12, 1)])
        self.answer = ([[2, 1], [1, 1]],
                       (2, 1),
                       (0, 0))

def suite():
    def is_test_class(member):
        name, value = member
        return inspect.isclass(value) and name.startswith("Test")

    classes = []
    members = inspect.getmembers(sys.modules[__name__])
    for name, value in members:
        if inspect.isclass(value) and name.startswith("Test"):
            classes.append(value)

    tests = map(unittest.TestLoader().loadTestsFromTestCase, classes)
    return unittest.TestSuite(tests)

if __name__ == "__main__":
    unittest.main()
