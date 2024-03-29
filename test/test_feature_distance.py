#!/usr/bin/env python
from __future__ import with_statement

import inspect
import os
import sys
import unittest

from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile

from segtools.feature_distance import DELIM, feature_distance, PRINT_FIELDS


def data2bed(lines):
    # Input lines should each be a tuple of (chrom, start, end, label)
    strings = []
    for line in lines:
        strings.append("chr%s\t%d\t%d\t%s" % line)
    return "\n".join(strings)

def data2gff(lines):
    # Input lines should each be a tuple of (chrom, start, end, label[, strand])
    strings = []
    for line in lines:
        source = "source"
        chr, start, end, label = line[0:4]
        start += 1
        chr = "chr%s" % chr

        fields = [chr, source, label, start, end]
        # Optionally add score and strand
        try:
            strand = line[4]
            fields.extend([".", strand])
        except IndexError:
            pass

        field_strs = [str(field) for field in fields]
        strings.append("\t".join(field_strs))

    return "\n".join(strings)

class OutputSaver(object):
    """Class to validate the stout behavior of a call

    Suggested usage:
      saver = OutputSaver()  # Start recording stdout
      <commands with stout recorded>
      saver.restore()  # Stop recording stdout
      log = saver.log()  # Retrieve stdout recording
    """
    def __init__(self):
        self._stdout = sys.stdout
        sys.stdout = self
        self._log = []
    def write(self, string):
        self._log.append(string)
    def writelines(self, lines):
        self._log.extend(lines)
    def restore(self):
        sys.stdout = self._stdout
    def log(self):
        return "".join(self._log)

class FeatureDistanceTestCases:
    class MainTester(unittest.TestCase):
        def init(self):
            pass

        def setUp(self):
            self.kwargs = {"verbose": False}
            self.init()

            # Set self.features from self.segments, self.features_list
            self.feature_file = None
            self.segment_file = None
            self._open_files = []
            
            self.output_path = mkdtemp()

            new_file = NamedTemporaryFile(mode="w",suffix=".bed")
            new_file.write(data2bed(self.segments))
            new_file.flush()
            self._open_files.append(new_file)
            self.segment_file = new_file.name

            type, data = self.features

            if type == "bed":
                new_file = NamedTemporaryFile(mode="w", suffix=".bed")
                new_file.write(data2bed(data))
            elif type == "gff":
                new_file = NamedTemporaryFile(mode="w",suffix=".gff")
                new_file.write(data2gff(data))
            elif type == "gtf":
                new_file = NamedTemporaryFile(mode="w",suffix=".gtf")
                new_file.write(data2gff(data))

            new_file.flush()
            self._open_files.append(new_file)
            self.feature_file = new_file.name

        def tearDown(self):
            for file in self._open_files:
                file.close()
            
            rmtree(self.output_path)

        def test(self):
            # Run function, logging stdout
            self._saver = OutputSaver()
            _ = feature_distance(self.segment_file, self.feature_file,
                                 self.output_path, **self.kwargs)
            self._saver.restore()
            log_lines = self._saver.log().strip("\n").split("\n")
            log_header = log_lines.pop(0).split(DELIM)

            # Check header
            header_answer = PRINT_FIELDS
            self.assertEqual(log_header, header_answer)

            # Check values
            observed_strs = [[line.split(DELIM)[-1]]
                            for line in log_lines]
            answer_strs = [[str(val) for val in line] for line in self.answer]

            self.assertEqual(observed_strs, answer_strs)

class TestSimple(FeatureDistanceTestCases.MainTester):
    def init(self):
        self.segments = [(1, 5, 6, 1),
                         (1, 10, 15, 1),
                         (1, 15, 20, 1)]
        self.features = ("bed", [(1, 6, 11, 1),
                                 (1, 11, 12, 1)])
        self.answer = [[1], [0], [4]]


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
