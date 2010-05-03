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
from numpy import array

from .bed import read_native as read_bed
from .common import check_clobber, DTYPE_SEGMENT_KEY, DTYPE_SEGMENT_START, DTYPE_SEGMENT_END, inverse_dict, maybe_gzip_open

PICKLED_EXT = "seg"
PICKLED_SUFFIX = os.extsep + PICKLED_EXT

class SegmentOverlapError(ValueError):
    pass

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

    def __init__(self, filename, verbose=True, pickled=None):
        """Returns a segmentation object derived from the given BED3+ file

        filename: path to the segmentation file (bed or pickled Segmentation)
        pickled: read filename as bed (False), pickled (True), or determine
          based upon the extension (None)
        """

        if not os.path.isfile(filename):
            raise IOError("Could not find Segmentation: %s" % filename)

        if verbose: print >>sys.stderr, "Loading segmentation:"

        if pickled is None:
            pickled = filename.endswith(PICKLED_SUFFIX)

        if pickled:
            # Read a pickled segmentation file in
            self._from_pickle(filename, verbose=verbose)
        else:
            # Parse a segmentation from a BED file
            self._from_bed(filename, verbose=verbose)


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

    def _from_pickle(self, filename, verbose=True):
        import cPickle

        if verbose: print >>sys.stderr, "  Unpickling Segmentation object...",
        with open(filename, 'rb') as ifp:
            self.__dict__ = cPickle.load(ifp).__dict__

        if verbose: print >>sys.stderr, "done"

    def pickle(self, namebase, verbose=True, clobber=False):
        """Pickle the segmentation into an output file"""
        import cPickle

        filename = namebase + PICKLED_SUFFIX

        check_clobber(filename, clobber)
        if verbose:
            print >>sys.stderr, ("Pickling Segmentation object to file: %s..."
                                 % filename),
        with open(filename, 'wb') as ofp:
            cPickle.dump(self, ofp, -1)

        if verbose: print >>sys.stderr, "done"

    def _from_bed(self, filename, verbose=True):
        """Parses a bedfile and sets some of the Segmentation fields to its data

        If the file is in BED3 format, labels will be None

        """
        metadata = self.get_bed_metadata(filename)

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

                assert datum.chromEnd >= datum.chromStart
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
        if unsorted_chroms:
            if verbose:
                print >>sys.stderr, ("  Segments were unsorted relative to the"
                                     " following chromosomes: %s" %
                                     ", ".join(unsorted_chroms))
                print >>sys.stderr, "  Sorting...",

            for chrom in unsorted_chroms:
                segments = chromosomes[chrom]
                segments.sort(order=['start'])

            if verbose: print >>sys.stderr, "done"

        if verbose:
            print >>sys.stderr, "  Checking for overlapping segments...",
        for chrom, segments in chromosomes.iteritems():
            if segments.shape[0] > 1:
                # Make sure there are no overlapping segments
                if (segments['end'][:-1] > segments['start'][1:]).any():
                    raise SegmentOverlapError("Found overlapping segments"
                                              " in chromosome: %s" % chrom)
        if verbose:
            print >>sys.stderr, "done"

        self.segtool, self.tracks = metadata
        self.name = filename
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

if __name__ == "__main__":
    pass

