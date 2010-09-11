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
from functools import partial
from numpy import array, int64, uint32

PICKLED_EXT = "pkl"

DTYPE_SEGMENT_START = int64
DTYPE_SEGMENT_END = int64
DTYPE_SEGMENT_KEY = uint32
DTYPE_STRAND = '|S1'

class Annotations(object):
    """Base class for representing annotations (BED/GFF/GTF files)

    chromosomes: a dict
      key: chromosome name
      val: segments, a numpy.ndarray, (each row is a (start, end, key) struct)
           sorted by start, end
           * These segments are not necessarily non-overlapping
    labels: a dict
      key: int ("label_key")  (a unique id; segment['key'])
      val: printable (str or int)  (what's in the actual file)
           This is the 4th column in a BED file and the 3rd column
           in GFF and GTF files.

    filename: the filename from which the annotations were loaded

    """
    class UnpickleError(Exception):
        pass

    class FilenameError(Exception):
        pass

    class FormatError(Exception):
        pass

    def __init__(self, filename, verbose=True):
        """Returns an Annotations object derived from the given file

        filename: path to a data file in one of the following formats:
          BED3+, GFF, GTF, or pickled Annotation. Format must be specified
          by the extension of the file (bed, gff, gtf, pkl), with all but
          pkl potentially gzipped (gz).
        """

        if not os.path.isfile(filename):
            raise IOError("Could not find file: %s" % filename)

        log("Loading %s:" % self.__class__, verbose)
        self._load(filename, verbose=verbose)

    @staticmethod
    def _get_file_format(filepath):
        """Determine the file format based upon the file extension"""
        filename = os.path.basename(filepath)
        root, ext = os.path.splitext(filename)

        if ext:
            # If gzip'd, process filename further
            if ext == ".gz":
                root, ext = os.path.splitext(root)

            if ext:
                return ext[1:].lower()  # remove dot

        raise Annotations.FilenameError("File missing appropriate extension:"
                                        " %s" % filepath)

    def _load(self, filename, verbose=True):
        format = self._get_file_format(filename)
        if format == PICKLED_EXT:
            # Read pickled object file
            self._from_pickle(filename, verbose=verbose)
        elif format in set(["bed", "gff", "gtf"]):
            self._from_file(filename, verbose=verbose)
        else:
            raise self.FormatError("Unrecognized file format: %s on file: %s"
                                   % (format, filename))

    def _iter_rows(self, filename, verbose=True):
        from .bed import read_native as read_bed
        from .gff import read_native as read_gff
        from .common import maybe_gzip_open

        format = self._get_file_format(filename)

        if format == "bed":
            reader = read_bed
        elif format == "gff":
            reader = read_gff
        elif format == "gtf":
            reader = partial(read_gff, gtf=True)
        else:
            raise self.FormatError("Cannot iterate segments in file format:"
                                   " %s" % format)

        with maybe_gzip_open(filename) as infile:
            for datum in reader(infile):
                row = {}
                d = datum.__dict__
                if format == "bed":
                    row['chrom'] = d['chrom']
                    row['start'] = d['chromStart']
                    row['end'] = d['chromEnd']
                    row['name'] = d.get('name', "")
                    row['strand'] = d.get('strand', ".")
                elif format == "gff" or format == "gtf":
                    row['chrom'] = d['seqname']
                    row['start'] = d['start']
                    row['end'] = d['end']
                    row['name'] = d.get('feature', "")
                    row['strand'] = d.get('strand', ".")

                if format == "gtf":
                    attrs = datum.attributes
                    row['gene_id'] = attrs['gene_id']
                    row['transcript_id'] = attrs['transcript_id']

                yield row

    def _from_pickle(self, filename, verbose=True):
        import cPickle

        log("  Unpickling %s object..." % self.__class__,
            verbose, end="")

        with open(filename, 'rb') as ifp:
            object = cPickle.load(ifp)
            if object.__class__ != self.__class__:
                raise self.UnpickleError("Error: Trying to load an indexed %s"
                                         " object as an indexed %s object!")
            self.__dict__ = object.__dict__

        log(" done", verbose)

    def _from_file(self, filename, verbose=True):
        """Parses a data file and sets object attributes

        Missing labels will be empty strings ("")

        """
        from .common import inverse_dict

        data = defaultdict(list)  # A dictionary-like object
        label_dict = {}
        last_segment_start = {}  # per chromosome
        unsorted_chroms = set()
        n_label_segments = {}
        n_label_bases = {}
        stranded = None
        for row in self._iter_rows(filename, verbose=verbose):
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            label = row['name']
            strand = row['strand']

            assert end >= start

            # Keep track of sorted chromosomes
            try:
                if start < last_segment_start[chrom]:
                    unsorted_chroms.add(chrom)
            except KeyError:
                pass

            # If any strands are specified, they all should be
            if strand == "+" or strand == "-":
                assert stranded is None or stranded
                stranded = True
            else:
                assert not stranded
                stranded = False

            try:  # Lookup label key
                label_key = label_dict[label]
            except KeyError:  # Map new label to key
                label_key = len(label_dict)
                label_dict[label] = label_key
                n_label_segments[label] = 0
                n_label_bases[label] = 0

            segment = [start, end, label_key]
            if stranded:
                segment.append(strand)

            data[chrom].append(tuple(segment))
            n_label_segments[label] += 1
            n_label_bases[label] += end - start
            last_segment_start[chrom] = start

        # Create reverse dict for field
        labels = inverse_dict(label_dict)

        # convert lists of tuples to array
        dtype = [('start', DTYPE_SEGMENT_START),
                 ('end', DTYPE_SEGMENT_END),
                 ('key', DTYPE_SEGMENT_KEY)]
        if stranded:
            dtype.append(('strand', DTYPE_STRAND))

        chromosomes = dict((chrom, array(segments, dtype=dtype))
                           for chrom, segments in data.iteritems())

        # Sort segments within each chromosome
        if unsorted_chroms:
            log("  Segments were unsorted relative to the following"
                " chromosomes: %s" % ", ".join(unsorted_chroms), verbose)
            log("  Sorting...", verbose, end="")

            for chrom in unsorted_chroms:
                segments = chromosomes[chrom]
                segments.sort(order=['start'])

            log(" done", verbose)

        self.filename = filename
        self.chromosomes = chromosomes
        self._labels = labels
        self._n_label_segments = n_label_segments
        self._n_label_bases = n_label_bases
        self._inv_labels = label_dict

    def pickle(self, verbose=True, clobber=False):
        """Pickle the annotations into an output file"""
        import cPickle
        from .common import check_clobber

        filename = os.extsep.join([self.filename, PICKLED_EXT])

        check_clobber(filename, clobber)
        log("Pickling %s object to file: %s..." % (self.__class__, filename),
            verbose, end="")
        with open(filename, 'wb') as ofp:
            cPickle.dump(self, ofp, -1)

        log(" done", verbose)

    def num_label_segments(self, label):
        return self._n_label_segments[str(label)]

    def num_segments(self):
        return sum(self._n_label_segments.values())

    def num_label_bases(self, label):
        return self._n_label_bases[str(label)]

    def num_bases(self):
        return sum(self._n_label_bases.values())

    def label_key(self, label):
        return self._inv_labels[str(label)]

    def label(self, label_key):
        return self._labels[label_key]

    @property
    def labels(self):
        return dict(self._labels)

    def set_labels(self, labels):
        missing = set(self._labels.keys()).difference(set(labels.keys()))

        if missing:
            raise ValueError("New labels do not cover old label keys: %r"
                             % missing)
        else:
            self._labels = labels

class Segmentation(Annotations):
    """Representation of a segmentation

    Segments must be non-overlapping.
    Strand information is not included.

    Defines additional attributes:
      tracks: a list of track names that were used to obtain the segmentation
      segtool: the tool that was used to obtain this segmentation (e.g. segway)
      * both will be empty if the file was not a BED file or did not contain
        this information
    """


    class SegmentOverlapError(ValueError):
        pass

    def _iter_rows(self, filename, verbose=True):
        """Override default row reading to ignore strand for Segmentations"""
        for row in super(self.__class__, self)._iter_rows(filename,
                                                          verbose=verbose):
            row['strand'] = '.';
            yield row;

    def _from_file(self, filename, verbose=True):
        """Wrap file loading to ensure no segments overlap

        This is preferred to wrapping __init__ because we don't need to check
        overlapping and don't want to load metadata if the Segmentation is
        being loaded from a pickle file instead of a BED/GFF file.
        """
        super(self.__class__, self)._from_file(filename, verbose=verbose)

        log("  Checking for overlapping segments...", verbose, end="")
        for chrom, segments in self.chromosomes.iteritems():
            # Make sure there are no overlapping segments
            if segments.shape[0] > 1 and \
                    (segments['end'][:-1] > segments['start'][1:]).any():
                raise self.SegmentOverlapError("Found overlapping segments"
                                               " in chromosome: %s" % chrom)

        log(" done", verbose)
        self.segtool, self.tracks = self.get_bed_metadata(filename)

    @staticmethod
    def get_bed_metadata(filename):
        from .common import maybe_gzip_open
        regexp = re.compile('description="(.*) segmentation of (.*)"')

        segtool = ""
        tracks = []

        with maybe_gzip_open(filename) as ifp:
            line = ifp.readline()

        matching = regexp.search(line)
        if matching:
            segtool = matching.group(1)
            tracks = matching.group(2).split(', ')

        return (segtool, tracks)


def log(message, verbose=True, end="\n", file=sys.stderr):
    """Wrapper for logging messages to stderr

    Similar to Python 3.0 print() syntax

    """
    if verbose:
        file.write(str(message))
        file.write(end)

if __name__ == "__main__":
    pass

