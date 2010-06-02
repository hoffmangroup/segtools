#!/usr/bin/env python
from __future__ import with_statement

import inspect
import sys
import unittest

from tempfile import NamedTemporaryFile

from segtools import Segmentation
from segtools.bed_compare import edit_distance

def beddata2bed(lines):
    # Input lines should each be a tuple of (chrom, start, end, label)
    strings = []
    for line in lines:
        strings.append("chr%s\t%d\t%d\t%s" % line)
    return "\n".join(strings)

class BedGenerator(unittest.TestCase):
    def init(self):
        pass

    def setUp(self):
        self.init()

        # Create fake bed files
        self._open_bedfiles = [NamedTemporaryFile(suffix=".bed"),
                               NamedTemporaryFile(suffix=".bed")]
        for file, data in zip(self._open_bedfiles, self.beddata):
            file.write(beddata2bed(data))
        for file in self._open_bedfiles:
            file.flush()
        self.bedfiles = [file.name for file in self._open_bedfiles]

    def tearDown(self):
        for file in self._open_bedfiles:
            file.close()

class EditDistanceTester(BedGenerator):
    def test(self):
        stats = edit_distance(self.bedfiles[0], self.bedfiles[1], verbose=False)
        self.assertValuesEqual(stats, self.answer)

    def assertValuesEqual(self, observed, expected):
        for val1, val2 in zip(observed, expected):
            self.assertEqual(val1, val2)

class TestPerfectOverlap(EditDistanceTester):
    def init(self):
        self.beddata = [[(1, 0, 50, 1),
                         (1, 50, 100, 2)],
                        [(1, 0, 50, 1),
                         (1, 50, 100, 2)]]
        self.answer = (100, 0, 0, 0)

class TestNoOverlap(EditDistanceTester):
    def init(self):
        self.beddata = [[(1, 10, 20, 1)],
                        [(1, 20, 35, 6)]]
        self.answer = (0, 0, 15, 10)

class TestSimpleDiff(EditDistanceTester):
    def init(self):
        self.beddata = [[(1, 11, 20, 1),
                         (1, 25, 35, 2)],
                        [(1, 15, 37, 2)]]
        self.answer = (10, 5, 7, 4)

class TestNotSegmentation(EditDistanceTester):
    def init(self):
        self.beddata = [[(1, 0, 10, 1),
                         (1, 5, 10, 1)],
                        [(1, 20, 35, 6)]]

    def test(self):
        self.assertRaises(Segmentation.SegmentOverlapError, edit_distance,
                          self.bedfiles[0], self.bedfiles[1], verbose=False)

def suite():
    classes = []
    members = inspect.getmembers(sys.modules[__name__])
    for name, value in members:
        if inspect.isclass(value) and name.startswith("Test"):
            classes.append(value)

    tests = map(unittest.TestLoader().loadTestsFromTestCase, classes)
    return unittest.TestSuite(tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
