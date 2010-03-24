#!/usr/bin/env python
from __future__ import division, with_statement

"""Segtools: tools for exploratory analysis of genomic segmentations

Copyright 2009: Orion Buske <stasis@uw.edu>

"""
__version__ = "$Revision$"

import os
import re
import sys

from collections import defaultdict
from pkg_resources import resource_listdir
from numpy import array

from .bed import read_native as read_bed
from .common import DTYPE_SEGMENT_KEY, DTYPE_SEGMENT_START, DTYPE_SEGMENT_END, inverse_dict, maybe_gzip_open

class Segmentation(object):
    """
    chromosomes: a dict
      key: chromosome name
      val: segments: numpy.ndarray, (each row is a (start, end, key) struct)
           sorted by start, end
           * These segments are not necessarily non-overlapping
    labels: a dict
      key: int ("label_key")  (a unique id; segment['end'])
      val: printable (str or int)  (what's in the actual BED file)

    tracks: a list of track names that were used to obtain the segmentation
    segtool: the tool that was used to obtain this segmentation (e.g. segway)
    name: the filename of the segmentation
    """

    def __init__(self, filename, verbose=True):
        """Returns a segmentation object derived from the given BED3+ file

        label_keys are integers corresponding to the order observed.
        """

        if not os.path.isfile(filename):
            raise IOError("Could not find Segmentation: %s" % filename)

        if verbose:
            print >>sys.stderr, "Loading segmentation...",

        # first get the tracks that were used for this segmentation
        self.segtool, self.tracks = self.get_bed_metadata(filename)

        self._load_bed(filename, verbose=verbose)

        if verbose:
    #         print >>sys.stderr, "\n\tMapped labels to integer keys:"
    #         for key, label in labels.iteritems():
    #             print >>sys.stderr, "\t\t\"%s\" -> %d" % (label, key)
            print >>sys.stderr, "done."

        self.name = filename

    @staticmethod
    def get_bed_metadata(filename):
        regexp = re.compile('description="(.*) segmentation of (.*)"')

        segtool = "Missing from BED file"
        tracks = ["Missing from BED file"]

        with maybe_gzip_open(filename) as ifp:
            line = ifp.readline()

        matching = regexp.search(line)
        if matching:
            segtool = matching.group(1)
            tracks = matching.group(2).split(', ')

        return (segtool, tracks)

    def _load_bed(self, filename, verbose=True):
        """Parses a bedfile and sets some of the Segmentation fields to its data

        If the file is in BED3 format, labels will be None

        """
        data = defaultdict(list)  # A dictionary-like object
        label_dict = {}
        last_segment_start = {}  # per chromosome
        unsorted_chroms = set()
        n_label_segments = {}
        n_label_bases = {}
        with maybe_gzip_open(filename) as infile:
            for datum in read_bed(infile):
                try:
                    if datum.chromStart < last_segment_start[datum.chrom]:
                        unsorted_chroms.add(datum.chrom)
                except KeyError:
                    pass

                try:
                    label = str(datum.name)
                except AttributeError:
                    # No name column, read as BED3 (no labels)
                    label = ""

                try:  # Lookup label key
                    label_key = label_dict[label]
                except KeyError:  # Map new label to key
                    label_key = len(label_dict)
                    label_dict[label] = label_key
                    n_label_segments[label] = 0
                    n_label_bases[label] = 0

                segment = (datum.chromStart, datum.chromEnd, label_key)
                data[datum.chrom].append(segment)
                last_segment_start[datum.chrom] = datum.chromStart

                n_label_segments[label] += 1
                n_label_bases[label] += segment[1] - segment[0]

        # Create reverse dict for field
        labels = inverse_dict(label_dict)

        # convert lists of tuples to array
        dtype = [('start', DTYPE_SEGMENT_START),
                 ('end', DTYPE_SEGMENT_END),
                 ('key', DTYPE_SEGMENT_KEY)]
        chromosomes = dict((chrom, array(segments, dtype=dtype))
                           for chrom, segments in data.iteritems())

        # Sort segments within each chromosome
        for chrom, segments in chromosomes.iteritems():
            if chrom in unsorted_chroms:
                if verbose:
                    print >>sys.stderr, "Segments were unsorted relative to \
    %s. Sorting..." % chrom,
                segments.sort(order=['start'])
                if verbose:
                    print >>sys.stderr, "done"

        self.chromosomes = chromosomes
        self._labels = labels
        self._n_label_segments = n_label_segments
        self._n_segments = sum(n_label_segments.values())
        self._n_label_bases = n_label_bases
        self._n_bases = sum(n_label_bases.values())
        self._inv_labels = label_dict
        return labels, chromosomes

    def num_label_segments(self, label):
        return self._n_label_segments[str(label)]

    def num_segments(self):
        return self._n_segments

    def num_label_bases(self, label):
        return self._n_label_bases[str(label)]

    def num_bases(self):
        return self._n_bases

    def label_key(self, label):
        return self._inv_labels[str(label)]

    def label(self, label_key):
        return self._labels[label_key]

    @property
    def labels(self):
        return dict(self._labels)

def test(verbose=False):
    import unittest

    # Gather a list of unittest modules
    filenames = resource_listdir(__name__, ".")
    regex = re.compile("^test_.*\.py$", re.IGNORECASE)
    module_filenames = filter(regex.search, filenames)
    def make_module_name(filename):
        return os.extsep.join([__name__, filename[:-3]])
    modulenames = map(make_module_name, module_filenames)
    print "Found test modules: %r" % modulenames
    map(__import__, modulenames)  # Import them all
    modules = [sys.modules[modulename] for modulename in modulenames]
    # Run the test suite for each
    suite = unittest.TestSuite([module.suite() for module in modules])
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    unittest.TextTestRunner(verbosity=verbosity).run(suite)


if __name__ == "__main__":
    pass

