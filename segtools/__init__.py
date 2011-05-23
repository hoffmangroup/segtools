#!/usr/bin/env python
from __future__ import division, with_statement

"""Segtools: tools for exploratory analysis of genomic segmentations

Copyright 2011: Michael Hoffman <mmh1@uw.edu>
Copyright 2009: Orion Buske <stasis@uw.edu>
"""
__version__ = "$Revision$"

import os
import re
import sys

from collections import defaultdict
from contextlib import contextmanager
from functools import partial
from numpy import array, int64, uint32
from time import time

from pkg_resources import resource_filename

try:
    PKG = __package__  # Python 2.6
except NameError:
    PKG = "segtools"

PKG_R = os.extsep.join([PKG, "R"])

EXT_GZ = "gz"
EXT_PICKLE = "pkl"
EXT_R = "R"

DTYPE_SEGMENT_START = int64
DTYPE_SEGMENT_END = int64
DTYPE_SEGMENT_KEY = uint32
DTYPE_STRAND = '|S1'

class Annotations(object):
    """Base class for representing annotations (BED/GFF/GTF files)

    chromosomes: a dict
      key: chromosome name
      val: segments, a numpy.ndarray,
           (each row is a (start, end, key, [strand, ...]) struct)
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

        log("Loading %s from file: %s" %
            (self.__class__.__name__, filename), verbose)
        start = time()
        self._load(filename, verbose=verbose)
        log("Loading finished in %.1f seconds." % (time() - start), verbose)

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
        if format == EXT_PICKLE:
            # Read pickled object file
            self._from_pickle(filename, verbose=verbose)
        elif format in set(["bed", "gff", "gtf"]):
            self._from_file(filename, verbose=verbose)
        else:
            raise self.FormatError("Unrecognized extension (%s) on file: %s"
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

        log("  Parsing lines from %s format" % format, verbose)

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
        from .common import maybe_gzip_open

        log("  Unpickling %s" % self.__class__.__name__, verbose)

        with maybe_gzip_open(filename, 'rb') as ifp:
            object = cPickle.load(ifp)
            if not issubclass(object.__class__, self.__class__):
                msg = ("Error: Cannot to load an indexed %s object"
                       " as an indexed %s object."
                       % (object.__class__.__name__,
                          self.__class__.__name__))
                raise self.UnpickleError(msg)

            self.__dict__ = object.__dict__

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
            log("  Sorting unsorted segments in the following"
                " chromosomes: %s" % ", ".join(unsorted_chroms), verbose)
            for chrom in unsorted_chroms:
                segments = chromosomes[chrom]
                segments.sort(order=['start'])

        self.filename = filename
        self.chromosomes = chromosomes
        self._labels = labels
        self._n_label_segments = n_label_segments
        self._n_label_bases = n_label_bases
        self._inv_labels = label_dict

    def pickle(self, filename=None, verbose=True, clobber=False):
        """Pickle the annotations into an output file"""
        import cPickle
        from .common import check_clobber, maybe_gzip_open

        if filename is None:
            filename = self.filename

        filename = os.extsep.join([filename, EXT_PICKLE, EXT_GZ])

        check_clobber(filename, clobber)
        log("Pickling %s object to file: %s"
            % (self.__class__.__name__, filename), verbose)
        with maybe_gzip_open(filename, 'wb') as ofp:
            cPickle.dump(self, ofp, -1)

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
        for row in super(self.__class__,
                         self)._iter_rows(filename, verbose=verbose):
            row['strand'] = '.';
            yield row;

    def _from_file(self, filename, verbose=True):
        """Wrap file loading to ensure no segments overlap

        This is preferred to wrapping __init__ because we don't need to check
        overlapping and don't want to load metadata if the Segmentation is
        being loaded from a pickle file instead of a BED/GFF file.
        """
        super(self.__class__, self)._from_file(filename, verbose=verbose)

        for chrom, segments in self.chromosomes.iteritems():
            # Make sure there are no overlapping segments
            if segments.shape[0] > 1 and \
                    (segments['end'][:-1] > segments['start'][1:]).any():
                raise self.SegmentOverlapError("Found overlapping segments"
                                               " in chromosome: %s" % chrom)

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

# XXX: should be replaced by use of logging module
# will allow percent-substitution only when verbose is on
def log(message, verbose=True, end="\n", file=sys.stderr):
    """Wrapper for logging messages to stderr

    Similar to Python 3.0 print() syntax

    """
    if verbose:
        file.write(str(message))
        file.write(end)

## Die with error message
def die(msg="Unexpected error"):
    log("\nFatal error: %s\n" % msg)
    sys.exit(1)


## Class to handle sourcing, plotting, and calling R functions
class RInterface(object):
    """
    Interface to R environment, providing methods to source files,
    call functions, and plot with R
    """
    _get_filename = partial(resource_filename, PKG_R)
    def __init__(self, files_to_source=[], verbose=True):
        self.started = False
        self._files = tuple(files_to_source)
        self.verbose = verbose
        self._r_console = lambda msg: sys.stderr.write(str(msg))

    def start(self, transcriptfile=None, verbose=None):
        if verbose is not None:
            self.verbose = verbose

        if self.started:
            return

        self.started = True

        # Start up R
        log("Opening connection to R environment.", self.verbose)
        from rpy2.robjects import r, rinterface
        self._r = r
        self._rinterface = rinterface
        self.RError = rinterface.RRuntimeError

        # Set up transcript
        self._transcript = transcriptfile
        if self._transcript:
            print >>self._transcript
            log("Saving transcript to file: %s" % self._transcript.name,
                self.verbose)

        # Source any R files
        for file in self._files:
            self.source(file)

    def source(self, filename):
        """Simplify importing R source in the package"""
        self.start()

        filename_full = self._get_filename(filename)
        if self._transcript:
            print >>self._transcript, "source(%r)" % filename_full

        try:
            self._r.source(filename_full)
        except self.RError:
            die("Failed to load R package: %s\n" % filename_full)

    @classmethod
    def arg_to_text(cls, arg):
        if isinstance(arg, bool):
            return repr(arg).upper() # True -> TRUE, False -> FALSE
        else:
            return repr(arg)

    def call(self, func, *args, **kwargs):
        """Safer way to call R functions (without importing all that junk)

        None values are substituted with empty string
        """
        self.start()

        from rpy2.robjects import numpy2ri
        # numpy2ri imported for side-effect

        # Make sure there are no None values in args or kwargs
        # rpy2 currently doesn't support passing Nones, so replace them with ""
        for arg_i, arg in enumerate(args):
            if arg is None:
                args[arg_i] = ""

        for key, val in kwargs.iteritems():
            if val is None:
                kwargs[key] = ""

        # Save call to transcript
        if self._transcript:
            args_text = [self.arg_to_text(arg) for arg in args]
            kwargs_text = ["%s = %s" % (key, self.arg_to_text(value))
                           for key, value in kwargs.iteritems()]
            all_args = ", ".join(args_text + kwargs_text)

            print >>self._transcript, "%s(%s)" % (func, all_args)

        return self._r[func](*args, **kwargs)

    def plot(self, func, *args, **kwargs):
        """Call R function but print timing information

        (kept around for backwards compatibility)
        """
        self.start()
        log("Plotting with R function: %r" % func, self.verbose)
        start = time()

        # Unless verbose, only print R output when there is an error
        r_log = []
        if self.verbose:
            r_console = self._r_console
        else:
            r_console = lambda msg: r_log.append(str(msg))

        self._rinterface.set_writeconsole(r_console)
        try:
            result = self.call(func, *args, **kwargs)
        except self.RError:
            die("Encoundered error within R:\n%s" % "\n  ".join(r_log))

        # Return console back to stderr
        self._rinterface.set_writeconsole(self._r_console)

        log("Plotting completed in %.1f seconds" % (time() - start),
            self.verbose)
        return result


@contextmanager
def open_transcript(dirpath, module, verified=False):
    filename = os.path.join(dirpath, os.extsep.join([module, EXT_R]))
    with open(filename, "w") as res:
        if not verified:
            print >>res, "## Experimental R transcript"
            print >>res, "## You may not be able to run the R code in this file exactly as written."
            print >>res

        yield res

if __name__ == "__main__":
    pass

