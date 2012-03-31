#!/usr/bin/env python
from __future__ import with_statement

import inspect
import sys
import unittest

from numpy import array
from tempfile import NamedTemporaryFile

from segtools import Annotation
from segtools.overlap import calc_overlap

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
        self.kwargs = {"verbose": False}
        self.init()

        # Set self.features from self.subject, self.query
        self.features = []
        self._open_files = []
        for type, data in [self.subject, self.query]:
            if type == "bed":
                new_file = NamedTemporaryFile(suffix=".bed")
                new_file.write(data2bed(data))
            elif type == "gff":
                new_file = NamedTemporaryFile(suffix=".gff")
                new_file.write(data2gff(data))
            elif type == "gtf":
                new_file = NamedTemporaryFile(suffix=".gtf")
                new_file.write(data2gff(data))

            new_file.flush()
            self._open_files.append(new_file)
            features = Annotation(new_file.name, verbose=False)

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

class TestSegmentPerfectOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 0, 50, 1),  # chrom, start, end, label
                                (1, 50, 100, 2)])
        self.query = ("bed", [(1, 0, 50, 1),
                              (1, 50, 100, 2)])
        self.answer = ([[1, 0], [0, 1]],  # subi vs qryj
                       (1, 1),  # total
                       (0, 0))  # none

class TestBasePerfectOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 0, 50, 1),
                                (1, 50, 100, 2)])
        self.query = ("bed", [(1, 0, 50, 1),
                              (1, 50, 100, 2)])
        self.answer = ([[50, 0], [0, 50]],
                       (50, 50),
                       (0, 0))

class TestSegmentNoOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 10, 20, 1)])
        self.query = ("bed", [(1, 20, 35, 6)])
        self.answer = ([0],
                       (1),
                       (1))

class TestBaseNoOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 10, 20, 1)])
        self.query = ("bed", [(1, 20, 35, 6)])
        self.answer = ([0],
                       (10),
                       (10))

class TestSegmentNoOverlapChrom(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 10, 20, 1)])
        self.query = ("bed", [(9, 20, 35, 6)])
        self.answer = ([0],
                       (1),
                       (1))

class TestBaseNoOverlapChrom(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 10, 20, 1)])
        self.query = ("bed", [(9, 20, 35, 6)])
        self.answer = ([0],
                       (10),
                       (10))

class TestSegmentSimpleOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 11, 20, 1),
                                (1, 25, 35, 2),
                                (1, 37, 40, 3)])
        self.query = ("gff", [(1, 15, 37, 2)])
        self.answer = ([[1], [1], [0]],
                       (1, 1, 1),
                       (0, 0, 1))

class TestBaseSimpleOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 11, 20, 1),
                                (1, 25, 35, 2),
                                (1, 37, 40, 3)])
        self.query = ("bed", [(1, 15, 37, 2)])
        self.answer = ([[5], [10], [0]],
                       (9, 10, 3),
                       (4, 0, 3))

class TestSegmentStackedFeatures(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 0, 10, 1)])
        self.query = ("gff", [(1, 2, 9, 2),
                              (1, 3, 12, 3),
                              (1, 4, 7, 4),
                              (1, 5, 6, 5)])
        self.answer = ([1, 1, 1, 1],
                       (1),
                       (0))

class TestBaseStackedFeatures(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 0, 10, 1)])
        self.query = ("gff", [(1, 2, 9, 2),
                              (1, 3, 12, 3),
                              (1, 4, 7, 4),
                              (1, 5, 6, 5)])
        self.answer = ([7, 7, 3, 1],
                       (10),
                       (2))

class TestSegmentOverlappingFeatures(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
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

class TestBaseOverlappingFeatures(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 0, 10, 1),
                                (1, 10, 20, 2),
                                (1, 20, 30, 1)])
        self.query = ("gff", [(1, 0, 30, 1),
                              (1, 2, 15, 1),
                              (1, 5, 20, 2),
                              (1, 10, 12, 1)])
        self.answer = ([[20, 5], [10, 10]],
                       (20, 10),
                       (0, 0))

class TestSegmentComplexOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 0, 4, 1),
                                (1, 10, 15, 2),
                                (1, 15, 21, 3),
                                (1, 26, 27, 1),
                                (1, 30, 31, 2),
                                (2, 14, 16, 2),
                                (3, 1, 100, 3)])
        self.query = ("bed", [(1, 0, 10, 1),
                              (1, 5, 15, 2),
                              (1, 15, 25, 1),
                              (1, 20, 30, 1),
                              (1, 20, 30, 2),
                              (1, 25, 30, 3),
                              (2, 5, 10, 1),
                              (2, 10, 15, 3)])
        self.answer = ([[2, 1, 1], [0, 1, 1], [1, 1, 0]],
                       (2, 3, 2),
                       (0, 1, 1))

class TestBaseComplexOverlap(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 0, 4, 1),
                                (1, 10, 15, 2),
                                (1, 15, 21, 3),
                                (1, 26, 27, 1),
                                (1, 30, 31, 2),
                                (2, 14, 16, 2),
                                (3, 1, 100, 3)])
        self.query = ("bed", [(1, 0, 10, 1),
                              (1, 5, 15, 2),
                              (1, 15, 25, 1),
                              (1, 20, 30, 1),
                              (1, 20, 30, 2),
                              (1, 25, 30, 3),
                              (2, 5, 10, 1),
                              (2, 10, 15, 3)])
        self.answer = ([[5, 1, 1], [0, 5, 1], [6, 1, 0]],
                       (5, 8, 105),
                       (0, 2, 99))

class TestSegmentComplexOverlapReversed(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "segments"
        self.subject = ("bed", [(1, 0, 10, 1),
                                (1, 5, 15, 2),
                                (1, 15, 25, 1),
                                (1, 20, 30, 1),
                                (1, 20, 30, 2),
                                (1, 25, 30, 3),
                                (2, 5, 10, 1),
                                (2, 10, 15, 3)])
        self.query = ("gff", [(1, 0, 4, 1),
                              (1, 10, 15, 2),
                              (1, 15, 21, 3),
                              (1, 26, 27, 1),
                              (1, 30, 31, 2),
                              (2, 14, 16, 2),
                              (3, 1, 100, 3)])
        self.answer = ([[2, 0, 2], [1, 1, 1], [1, 1, 0]],
                       (4, 2, 2),
                       (1, 0, 0))


class TestBaseComplexOverlapReversed(OverlapTester):
    def init(self):
        self.kwargs["mode"] = "bases"
        self.subject = ("bed", [(1, 0, 10, 1),
                                (1, 5, 15, 2),
                                (1, 15, 25, 1),
                                (1, 20, 30, 1),
                                (1, 20, 30, 2),
                                (1, 25, 30, 3),
                                (2, 5, 10, 1),
                                (2, 10, 15, 3)])
        self.query = ("gff", [(1, 0, 4, 1),
                              (1, 10, 15, 2),
                              (1, 15, 21, 3),
                              (1, 26, 27, 1),
                              (1, 30, 31, 2),
                              (2, 14, 16, 2),
                              (3, 1, 100, 3)])
        self.answer = ([[5, 0, 7], [1, 5, 1], [1, 1, 0]],
                       (35, 20, 10),
                       (23, 13, 8))

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
