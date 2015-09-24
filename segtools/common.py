from __future__ import division, with_statement

"""
Assorted utility functions and classes common or useful to most of Segtools.


Author: Orion Buske <stasis@uw.edu>
"""

# XXX: the things in here should be moved to __init__; test

import os

from contextlib import closing, contextmanager
from functools import partial
from gzip import open as _gzip_open
from warnings import warn
from time import time

from numpy import array, empty, iinfo, ndarray, searchsorted
from tabdelim import DictReader, DictWriter, ListReader, ListWriter

from . import log, EXT_GZ, PKG, die

PKG_RESOURCE = os.extsep.join([PKG, "resources"])

EXT_PDF = "pdf"
EXT_PNG = "png"
EXT_TAB = "tab"
EXT_DOT = "dot"
EXT_HTML = "html"
EXT_DIV = "div"
EXT_SLIDE = os.extsep.join(["slide", EXT_PNG])
EXT_THUMB = os.extsep.join(["thumb", EXT_PNG])  # for thumbnails
EXT_SUMMARY = os.extsep.join(["summary", EXT_TAB])

IMG_EXTS = [EXT_PNG, EXT_PDF, EXT_SLIDE, EXT_THUMB]
NICE_EXTS = dict(tab=EXT_TAB, pdf=EXT_PDF, png=EXT_PNG, html=EXT_HTML,
                 div=EXT_DIV, slide=EXT_SLIDE, thumb=EXT_THUMB,
                 summary=EXT_SUMMARY, dot=EXT_DOT)

SUFFIX_GZ = os.extsep + EXT_GZ

LABEL_ALL = "all"

THUMB_SIZE = 100

class OutputMasker(object):
    """Class to mask the output of a stream.

    Suggested usage:
      sys.stout = OutputMasker(sys.stdout)  # Start masking
      <commands with stout masked>  # Masked commands
      sys.stdout = sys.stdout.restore()  # Stop masking
    """
    def __init__(self, stream=None):
        self._stream = stream
    def write(self, string):
        pass  # Mask output
    def writelines(self, lines):
        pass  # Mask output
    def restore(self):
        return self._stream


class FeatureScanner(object):
    def __init__(self, features):
        self._current = None
        self._nearest = None
        if isinstance(features, ndarray):
            self._iter = iter(features)
            try:
                self._current = self._iter.next()
            except StopIteration:
                pass
            else:
                self._nearest = self._current
        else:
            warn("Feature list not of understood type. Treated as empty")

        self._has_next = (self._current is not None)
        self._overlapping = []

        self._prev_best = None
        self._next_best = None
        self._prev_start = None

    @staticmethod
    def distance(feature1, feature2):
        """Return the distance between two features (zero if overlap)"""
        if feature1['end'] <= feature2['start']:
            return feature2['start'] - feature1['end'] + 1
        elif feature2['end'] <= feature1['start']:
            return feature1['start'] - feature2['end'] + 1
        else:  # Features overlap
            return 0

    @staticmethod
    def are_overlapping(feature1, feature2):
        """Return the whether or not two features overlap"""
        return feature1['start'] < feature2['end'] and \
                feature2['start'] < feature1['end']

    def _seek(self, segment):
        # Seek internal records up to segment, updated self._overlapping
        # and self._nearest

        # Clear out any non-overlapping features from self._overlapping
        if self._overlapping:
            i = 0
            while i < len(self._overlapping):
                if self.are_overlapping(segment, self._overlapping[i]):
                    i += 1
                else:
                    del self._overlapping[i]

        # If any are overlapping, nearest is 0 distance
        if self._overlapping:
            nearest_dist = 0
            self._nearest = self._overlapping[0]
        elif self._nearest is not None:
            # Closest feature will start as the closest from the last iteration
            nearest_dist = self.distance(self._nearest, segment)
        else:
            assert not self._has_next

        # Scan through features until we get past segment
        if self._has_next:
            segment_start = segment['start']
            segment_end = segment['end']
            while self._current['start'] < segment_end:
                # dist will be negative if this overlaps
                dist = segment_start - self._current['end'] + 1
                if dist < 0:
                    self._overlapping.append(self._current)

                if dist < nearest_dist:
                    self._nearest = self._current
                    nearest_dist = dist

                try:
                    self._current = self._iter.next()
                except StopIteration:
                    self._has_next = False
                    break

            # Check distance of first feature past segment as well
            dist = self._current['start'] - segment_end + 1
            if dist < nearest_dist:
                self._nearest = self._current
                nearest_dist = dist

    def nearest(self, segment):
        """Return a list of the nearest features for the next segment

        This will contain a single feature if no features overlap the given
        segment, else it will be a list of features that overlap.

        If no features were found at all, this will be None.

        Should be called with in-order segments

        """
        self._seek(segment)
        if self._overlapping:
            return self._overlapping
        elif self._nearest is None:
            return None
        else:
            return self._nearest


## UTILITY FUNCTIONS

def inverse_dict(d):
    """Given a dict, returns the inverse of the dict (val -> key)"""
    res = {}
    for k, v in d.iteritems():
        assert v not in res
        res[v] = k
    return res

def make_filename(dirpath, basename, ext):
    return os.path.join(dirpath, os.extsep.join([basename, ext]))

make_tabfilename = partial(make_filename, ext=EXT_TAB)
make_htmlfilename = partial(make_filename, ext=EXT_HTML)
make_divfilename = partial(make_filename, ext=EXT_DIV)
make_pngfilename = partial(make_filename, ext=EXT_PNG)
make_pdffilename = partial(make_filename, ext=EXT_PDF)
make_thumbfilename = partial(make_filename, ext=EXT_THUMB)
make_slidefilename = partial(make_filename, ext=EXT_SLIDE)
make_summaryfilename = partial(make_filename, ext=EXT_SUMMARY)
make_dotfilename = partial(make_filename, ext=EXT_DOT)

def make_namebase_summary(namebase):
    return os.extsep.join([namebase, "summary"])

def make_id(modulename, dirpath):
    return "_".join([modulename, os.path.basename(dirpath)])

def check_clobber(filename, clobber):
    if (not clobber) and os.path.isfile(filename):
        raise IOError("Output file: %s already exists."
                      " Use --clobber to overwrite." % filename)

def gzip_open(*args, **kwargs):
    return closing(_gzip_open(*args, **kwargs))

def maybe_gzip_open(filename, *args, **kwargs):
    if filename.endswith(SUFFIX_GZ):
        return gzip_open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)

def fill_array(scalar, shape, dtype=None, *args, **kwargs):
    if dtype is None:
        dtype = array(scalar).dtype

    res = empty(shape, dtype, *args, **kwargs)
    res.fill(scalar)
    return res

@contextmanager
def tab_saver(dirpath, basename, fieldnames=None, header=True,
              clobber=False, metadata=None, verbose=True):
    """Save data to tab file

    If fieldnames are specified, a DictWriter is yielded
    If fieldnames are not, a ListWriter is yielded instead

    metadata: a dict describing the data to include in the comment line.
      Comment line will start with '# ' and then will be a space-delimited
      list of <field>=<value> pairs, one for each element of the dict.
    """
    outfilename = make_tabfilename(dirpath, basename)
    log("Saving tab file: %s" % outfilename, verbose)
    start = time()
    check_clobber(outfilename, clobber)
    with open(outfilename, "w") as outfile:
        if metadata is not None:
            assert isinstance(metadata, dict)
            metadata_strs = ["%s=%s" % pair for pair in metadata.iteritems()]
            print >>outfile, "# %s" % " ".join(metadata_strs)

        if fieldnames:
            yield DictWriter(outfile, fieldnames, header=header,
                             extrasaction="ignore")
        else:
            yield ListWriter(outfile)

    log("Saved tab file in %.1f seconds" % (time() - start), verbose)

@contextmanager
def tab_reader(dirpath, basename, verbose=True, fieldnames=False):
    """Reads data from a tab file

    Yields a tuple (Reader, metadata)
    If fieldnames is True, a DictReader is yielded
    If fieldnames is False, a ListReader is yielded instead

    metadata: a dict describing the data included in the comment line.
    """

    infilename = make_tabfilename(dirpath, basename)
    log("Reading tab file: %s" % infilename, verbose)
    start = time()
    if not os.path.isfile(infilename):
        raise IOError("Unable to find tab file: %s" % infilename)

    with open(infilename, "rU") as infile:
        infile_start = infile.tell()
        comments = infile.readline().strip().split()
        metadata = {}
        if comments and comments[0] == "#":
            # Found comment line
            for comment in comments[1:]:
                name, value = comment.split("=")
                metadata[name] = value
        else:
            infile.seek(infile_start)

        if fieldnames:
            yield DictReader(infile), metadata
        else:
            yield ListReader(infile), metadata

    log("Tab file read in %.1f seconds" % (time() - start), verbose)

def map_segment_label(segments, range):
    """Flatten segments to a segment label per base vector

    Returns a tuple of:
      numpy.array of the segment key at every base
      sentinel value assigned to non-overlapped bases

    The given range must either be a tuple of (start, end), or
      have 'start' and 'end' attributes, and the returned array
      will represent the segment keys found between start and end.

    Thus, if start is 1000, and a segment with key 4 starts at 1000, the
    value of keys[0] is 4.

    """
    if isinstance(range, tuple):
        map_start, map_end = range
    else:
        map_start = range.start
        map_end = range.end

    map_size = map_end - map_start
    assert map_size >= 0

    # Choose sentinel value as maximum value supported by datatype
    segments_dtype = segments['key'].dtype
    sentinel = iinfo(segments_dtype).max
    assert sentinel != segments['key'].max()  # not already used

    res = fill_array(sentinel, (map_size,), segments_dtype)
    if map_size == 0:
        return res

    # will overwrite in overlapping case
    for segment in segments:
        start = segment['start']
        end = segment['end']
        key = segment['key']
        if start < map_end and end > map_start:
            start_i = max(start - map_start, 0)
            end_i = min(end - map_start, map_size)
            res[start_i:end_i] = key

    return res, sentinel

def iter_segments_continuous(chromosome, segmentation, column=None,
                             verbose=True):
    """Yields segment and the continuous corresponding to it, for each segment
    in the chromosome inside of a supercontig
    column: only report continuous from the specified column (or None for all)
    """
    chrom = chromosome.name
    try:
        segments = segmentation.chromosomes[chrom]
    except KeyError:
        raise StopIteration

    supercontig_iter = chromosome.itercontinuous()
    supercontig = None
    supercontig_last_start = 0

    for segment in segments:
        start = segment['start']
        end = segment['end']

        while supercontig is None or start >= supercontig.end:
            try:
                # Raises StopIteration if out of supercontigs
                supercontig, continuous = supercontig_iter.next()
            except StopIteration:
                raise

            # Enforce increasing supercontig indices
            assert supercontig.start >= supercontig_last_start
            supercontig_last_start = supercontig.start

            # If only pulling specific column, load entire supercontig
            # of data into memory
            if start < supercontig.end and column is not None:
                continuous = continuous[:, column]

        if end <= supercontig.start:
            # Jump forward to next segment
            continue

        sc_start = supercontig.project(start)
        sc_end = supercontig.project(end)
        yield segment, continuous[sc_start:sc_end]


## Yields supercontig and the subset of segments which overlap it.
def iter_supercontig_segments(chromosome, segmentation, verbose=True):
    try:
        segments = segmentation.chromosomes[chromosome.name]
    except KeyError:
        raise StopIteration

    for supercontig in chromosome:
        # Since it's a segmentation, there are no overlapping segments.
        # Thus, sorting by the start also sorts by the end, and we can
        # binary search on both to get the bounds of the overlapping segments.
        # Overlapping segments are guaranteed to be a contiguous set of rows
        # in the segment array for this reason.
        start = searchsorted(segments['end'], supercontig.start, side="right")
        end = searchsorted(segments['start'], supercontig.end, side="left")
        if end > start:
            yield supercontig, segments[start:end]

## Ensures output directory exists and has appropriate permissions
def setup_directory(dirname):
    if not os.path.isdir(dirname):
        try:
            os.makedirs(dirname)
        except IOError:
            die("Error: Could not create output directory: %s" % (dirname))
    else:
        # Require write and execute permissions on existing dir
        if not os.access(dirname, os.W_OK | os.X_OK):
            die("Error: Output directory is not writeable and executable.")

## Given labels and mnemonics, returns an ordered list of label_keys
##   and a new labels dict mapping label_keys to label strings
## If no mnemonics are specified, returns the passed labels and
##   a label_key ordering
def get_ordered_labels(labels, mnemonics=[]):
    if mnemonics is not None and len(mnemonics) > 0:
        # Create key lookup dictionary
        key_lookup = {}  # old_label -> label_key
        for label_key, label in labels.iteritems():
            assert(label not in key_lookup)  # Enforce 1-to-1
            key_lookup[label] = label_key

        labels = {}
        ordered_keys = []
        for old_label, new_label in mnemonics:
            label_key = key_lookup[old_label]
            ordered_keys.append(label_key)
            labels[label_key] = new_label
    else:
        # Don't change labels, but specify ordering
        partial_int_labels = {}
        for key, label in labels.iteritems():
            try:
                partial_int_labels[int(label)] = key
            except ValueError:
                partial_int_labels[label] = key
        ordered_keys = list(partial_int_labels[key]
                            for key in sorted(partial_int_labels.keys()))

    return ordered_keys, labels

## Maps the provided labels to mnemonics (or some other mnemonic file
##   field specified ready to be fed into R.
## Returns a numpy.array of strings with a row of [oldlabel, newlabel] for
## every old_label, and their ordering specifies the desired display order
## Labels should be a dict of label_keys -> label strings
def map_mnemonics(labels, mnemonicfilename, field="new"):
    if not mnemonicfilename:
        return []

    label_order, label_mnemonics = load_mnemonics(mnemonicfilename)
    str_labels = labels.values()
    mnemonics = []
    # Add mapping for labels in mnemonic file
    for old_label in label_order:
        #try:
        new_label = label_mnemonics[old_label][field]
        # XXX: need KeyError for relabel.relabel()
        #except KeyError:
        #die("Mnemonic file missing expected column: %s" % field)

        if old_label in str_labels:
            mnemonics.append([old_label, new_label])

    # Add mapping for labels not in mnemonic file
    # Use ordering of label mapping
    for label_key in sorted(labels.keys()):
        label = labels[label_key]  # Get actual label (string)
        if label not in label_order:  # Then map to same name
            mnemonics.append([label, label])

    return array(mnemonics)

def make_comment_ignoring_dictreader(iterable, *args, **kwargs):
    return DictReader((item for item in iterable if not item.startswith("#")),
                      *args, **kwargs)

## Loads segmentation label descriptions and mnemonics
##   from a tab-delimited file with a header row
def load_mnemonics(filename):
    """
    Input file should have a tab-delimited row for each label, of the form:
               old    new    description      color
      e.g.     4    IS   initiation (strong)  255,0,0
    Returns a tuple of (label_order, label_mnemonics)

    label_order: an ordered list of old labels,
      corresponding to the preferred display order in plots

    label_mnemonics: dict
      key: a string corresponding to an "old" label
      value: a dict with fields including potentially "new", "description",
        and "color" where description is a several-word description of the
        label, new is a few-character identifier, and color is a UCSC Genome
        Browser integer RGB color triplet
    """
    if filename is None:
        return []
    elif not os.path.isfile(filename):
        die("Could not find mnemonic file: %s" % filename)

    label_order = []
    label_mnemonics = {}
    with open(filename, "rU") as ifp:
        reader = make_comment_ignoring_dictreader(ifp)
        for row in reader:
            try:
                label_index = row["old"]
            except KeyError:
                die("Mnemonic file missing required column: 'old'")

            label_order.append(label_index)
            label_mnemonics[label_index] = row

    return (label_order, label_mnemonics)
