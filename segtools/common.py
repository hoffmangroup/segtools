from __future__ import division, with_statement

__version__ = "$Revision$"

"""
Assorted utility functions and classes common or useful to most of segtools.


Author: Orion Buske <stasis@uw.edu>
"""

import os
import re
import sys

from collections import defaultdict
from contextlib import closing, contextmanager
from functools import partial
from genomedata import Genome
from gzip import open as _gzip_open
from numpy import array, concatenate, empty, iinfo, logical_and, uint16, uint32, int64
from operator import itemgetter
from pkg_resources import resource_filename
from tabdelim import DictReader, DictWriter, ListReader, ListWriter
from rpy2.robjects import r, rinterface, numpy2ri
# numpy2ri imported for side-effect

try:
    PKG = __package__  # Python 2.6
except NameError:
    PKG = "segtools"

PKG_R = os.extsep.join([PKG, "R"])
PKG_RESOURCE = os.extsep.join([PKG, "resources"])

from .gff import read_native as read_gff

EXT_GZ = "gz"
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

DTYPE_SEGMENT_START = int64
DTYPE_SEGMENT_END = int64
DTYPE_SEGMENT_KEY = uint32
DTYPE_SOURCE_KEY = uint16
DTYPE_STRAND = '|S1'
MAXLEN_ID_DTYPE = 20
DTYPE_GENE_ID = '|S%d' % MAXLEN_ID_DTYPE
DTYPE_TRANSCRIPT_ID = '|S%d' % MAXLEN_ID_DTYPE

#FEATURE_FIELDS = ["chrom", "source", "name", "start", "end", "score", "strand"]

r_filename = partial(resource_filename, PKG_R)


## Wrapper for gff/gtf feature data
class Features(object):
    """
    sequences: a dict
      key: chromosome/scaffold name
      val: segments: numpy.ndarray
           each row is a struct of: start, end, key, source_key
           sorted by start, end
    features: a dict
      key: int ("feature_key")  (a unique id)
      val: printable (str or int)  (what's in the actual GFF/GTF file)
    sources: a dict
      key: int ("source_key")  (a unique id)
      val: printable (str or int)  (what's in the actual GFF/GTF file)

    tracks: a list of track names that were used to obtain the segmentation
    segtool: the tool that was used to obtain this segmentation (e.g. segway)
    """
    def __init__(self, sequences, features, sources):
        self.sequences = sequences
        self.features = features
        self.sources = sources
        # Allow Segmentation-like access
        self.chromosomes = self.sequences
        self.labels = self.features

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

class ProgressBar(object):
    def __init__(self, total=67, width=67, chr_done="#", chr_undone="-",
                 format_str="Progress: [%(done)s%(undone)s]",
                 out=sys.stdout):
        """Create a progress bar for num_items that is width characters wide

        total: the number of items in the task (calls to next before 100%)
        width: the width of the bar itself, not including labels
        chr_done: character for displaying completed work
        chr_undone: character for displaying uncompleted work
        format_str: format string for displaying progress bar. Must contain
          a reference to %(done)s and %(undone)s, which will be substituted
        out: an object that supports calls to write() and flush()
        """
        assert int(width) > 0
        assert int(total) > 0

        self._width = int(width)
        self._n = int(total)
        self._quantum = float(self._n / self._width)
        self._chr_done = str(chr_done)
        self._chr_undone = str(chr_undone)
        self._format_str = str(format_str)
        self._out = out

        self._i = 0
        self._progress = 0
        self.refresh()

    def next(self):
        """Advance to the next item in the task

        Might or might not refresh the progress bar
        """
        self._i += 1
        if self._i > self._n:
            raise StopIteration("End of progress bar reached.")

        progress = int(self._i / self._quantum)
        if progress != self._progress:
            self._progress = progress
            self.refresh()

    def refresh(self):
        """Refresh the progress bar display"""
        fields = {"done": self._chr_done * self._progress,
                  "undone": self._chr_undone * (self._width - self._progress)}
        self._out.write("\r%s" % (self._format_str % fields))
        self._out.flush()

    def end(self):
        """Complete the job progress, regardless of current state"""
        self._progress = self._width
        self.refresh()
        self._out.write("\n")
        self._out.flush()


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

## Die with error message
def die(msg="Unexpected error"):
    print >> sys.stderr, "\nERROR: %s\n" % msg
    sys.exit(1)

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
        raise IOError("Output file: %s already exists!"
                      " Use --clobber to overwrite!" % filename)

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


## XXX: No known usage
@contextmanager
def none_contextmanager():
    yield None

def r_source(filename):
    """
    Simplify importing R source in the package
    """
    try:
        r.source(r_filename(filename))
    except rinterface.RRuntimeError:
        die("Failed to load R package: %s\n" % filename)

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
    if verbose:
        print >>sys.stderr, "Saving tab file...",

    outfilename = make_tabfilename(dirpath, basename)
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

    if verbose:
        print >>sys.stderr, "done"

@contextmanager
def tab_reader(dirpath, basename, verbose=True, fieldnames=False):
    """Reads data from a tab file

    Yields a tuple (Reader, metadata)
    If fieldnames is True, a DictReader is yielded
    If fieldnames is False, a ListReader is yielded instead

    metadata: a dict describing the data included in the comment line.
    """
    if verbose:
        print >>sys.stderr, "Reading tab file...",

    infilename = make_tabfilename(dirpath, basename)
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

    if verbose:
        print >>sys.stderr, "done"

@contextmanager
def image_saver(dirpath, basename, clobber=False, verbose=True,
                downsample=False, **kwargs):
    """
    Generator to save an R plot to both png and pdf with only one plot
    Yields to caller to plot, then saves images
    """
    if verbose:
        print >>sys.stderr, "Creating images...",

    png_filename = make_pngfilename(dirpath, basename)
    check_clobber(png_filename, clobber)
    try:
        r.png(png_filename, **kwargs)
        png_device = r["dev.cur"]()
        r["dev.control"]("enable")
    except rinterface.RRuntimeError:
        die('Image creation failed.\nIf unable to start PNG device, try'
            ' setting (export/setenv) variable DISPLAY to ":1" from the'
            ' (bash/c) shell before re-running validation.')

    yield  # Wait for plot

    # Use R routine to create all the other images (pdf, slide, etc)
    try:
        r["dev.print.images"](basename=basename, dirname=dirpath,
                              downsample=downsample, clobber=clobber, **kwargs)
        r["dev.off"](png_device)
    except rinterface.RRuntimeError:
        print >> sys.stderr, 'ERROR: Images might not have been saved!'

    if verbose:
        print >>sys.stderr, "done"


## Maps label_keys over segments
def map_segments(segments, labels, chrom_size):
    """
    converts a segment mapping into label_keys at every nucleotide position
        in the chromosome

    e.g.  segment labeled 0 at position [0,4) followed by
          segment labeled 1 at position [4,7) gets converted into:
          0000111
    """
    segments_dtype = segments['key'].dtype  # MA: the data type of segments
    segments_dtype_max = iinfo(segments_dtype).max  # sentinel value
    # MA: the maximum value supported by the type
    assert segments_dtype_max not in labels.keys()

    res = fill_array(segments_dtype_max, (chrom_size,), segments_dtype)

    # will overwrite in overlapping case
    for start, end, key in segments:
        res[start:end] = key

    return res

def map_segment_label(segments, range):
    """Flatten segments to a segment label per base vector

    Returns a numpy.array of the segment key at every base

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

    # Choose sentinal value as maximum value supported by datatype
    segments_dtype = segments['key'].dtype
    sentinal = iinfo(segments_dtype).max
    assert sentinal != segments['key'].max()  # not already used

    res = fill_array(sentinal, (map_size,), segments_dtype)
    if map_size == 0:
        return res

    # will overwrite in overlapping case
    for start, end, key in segments:
        if start < map_end and end > map_start:
            start_i = max(start - map_start, 0)
            end_i = min(end - map_start, map_size)
            res[start_i:end_i] = key

    return res

## Yields segment and the continuous corresponding to it, for each segment
##   in the chromosome inside of a supercontig
def iter_segments_continuous(chromosome, segmentation, verbose=True):
    chrom = chromosome.name
    try:
        segments = segmentation.chromosomes[chrom]
    except KeyError:
        raise StopIteration

    supercontig_iter = chromosome.itercontinuous()
    supercontig = None
    supercontig_last_start = 0
    nsegments = len(segments)

    if verbose:
        format_str = "".join([chrom, ":\t[%(done)s%(undone)s]"])
        progress = ProgressBar(nsegments, out=sys.stderr,
                               format_str=format_str)

    for segment in segments:
        start = segment['start']
        end = segment['end']

        while supercontig is None or start >= supercontig.end:
            try:
                # Raise StopIteration if out of supercontigs
                supercontig, continuous = supercontig_iter.next()
            except StopIteration:
                if verbose:
                    progress.end()
                raise

            # Enforce increasing supercontig indices
            assert supercontig.start >= supercontig_last_start
            supercontig_last_start = supercontig.start

        if verbose:
            progress.next()

        if end <= supercontig.start:
            continue  # Get next segment

        try:
            sc_start = supercontig.project(max(start, supercontig.start))
            sc_end = supercontig.project(min(end, supercontig.end))
            yield segment, continuous[sc_start:sc_end]
        except:
            for k, v in locals():
                print >>sys.stderr, "%r: %r" % (k, v)
            raise

    if verbose:
        progress.end()

## Yields supercontig and the subset of segments which overlap it.
def iter_supercontig_segments(chromosome, segmentation, verbose=True):
    try:
        segments = segmentation.chromosomes[chromosome.name]
    except KeyError:
        raise StopIteration

    for supercontig in chromosome:
        if verbose:
            print >>sys.stderr, "\t\t%s" % supercontig

        rows = logical_and(segments['start'] < supercontig.end,
                           segments['end'] > supercontig.start)
        cur_segments = segments[rows]
        if cur_segments.shape[0] > 0:
            yield supercontig, segments[rows]

## Ensures output directory exists and has appropriate permissions
def setup_directory(dirname):
    if not os.path.isdir(dirname):
        try:
            os.mkdir(dirname)
        except IOError:
            die("Error: Could not create output directory: %s" % (dirname))
    else:
        # Require write and execute permissions on existing dir
        if not os.access(dirname, os.W_OK | os.X_OK):
            die("Error: Output directory does not have appropriate"
                " permissions!")

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

    label_order, label_data = load_mnemonics(mnemonicfilename)
    str_labels = labels.values()
    mnemonics = []
    # Add mapping for labels in mnemonic file
    for old_label in label_order:
        try:
            new_label = label_data[old_label][field]
        except KeyError:
            die("Mnemonic file missing expected column: %s" % field)

        if old_label in str_labels:
            mnemonics.append([old_label, new_label])

    # Add mapping for labels not in mnemonic file
    # Use ordering of label mapping
    for label_key in sorted(labels.keys()):
        label = labels[label_key]  # Get actual label (string)
        if label not in label_order:  # Then map to same name
            mnemonics.append([label, label])

    return array(mnemonics)

## Loads segmentation label descriptions and mnemonics
##   from a tab-delimited file with a header row
def load_mnemonics(filename):
    """
    Input file should have a tab-delimited row for each label, of the form:
               old    new    description
      e.g.     4    IS   initiation (strong)
    Returns a tuple of (label_order, label_data)

    label_order: an ordered list of old labels,
      corresponding to the preferred display order in plots

    label_mnemonics: dict
      key: a string corresponding to an "old" label
      value: a dict with fields including "new" and "description",
        where description is a several-word description of the label
        and new is a few-character identifier
    """
    if filename is None:
        return []
    elif not os.path.isfile(filename):
        die("Could not find mnemonic file: %s" % filename)

    label_order = []
    label_data = {}
    with open(filename, "rU") as ifp:
        reader = DictReader(ifp)
        for row in reader:
            try:
                label_index = row["old"]
            except KeyError:
                die("Mnemonic file missing required column: 'old'")

            label_order.append(label_index)
            label_data[label_index] = row

    return (label_order, label_data)

## Parses gff file for features
def load_gff_data(gff_filename, sort=True):
    '''
    Expects data in GFF format (1-indexed locations):
    CHROM<tab>source<tab>feature<tab>START<tab>END<tab>score<tab>STRAND
    chrom, start, end, and strand are required
    strand can be '+'/'-' or '.', but not both in the same file

    File may be gzip'd, but if so, must have .gz as the extension

    Returns gffdata
    gffdata: a dict
      key: string chromosome name (e.g. "chr13")
      val: a list of dicts (ordered by the start index, ascending)
        key: "start", "end", "strand", "source"
        val: the associated data item, if it exists
             string for strand
             int (zero-indexed) for start and end (exclusive)
    '''
    data = defaultdict(list)
    stranded = None
    with maybe_gzip_open(gff_filename) as infile:
        for line in infile:
            # Ignore comments
            if line.startswith("#"): continue

            # Parse tokens from GFF line
            try:
                fields = {}
                tokens = line.strip().split("\t")

                chrom = tokens[0]
                fields["source"] = tokens[1]
                fields["name"] = tokens[2]
                fields["start"] = int(tokens[3]) - 1  # Make zero-indexed
                # Make zero-indexed and exclusive (these cancel out)
                fields["end"] = int(tokens[4])

                try:
                    strand = tokens[6]
                except IndexError:
                    strand = "."

                if strand == "+" or strand == "-":
                    assert stranded or stranded is None
                    stranded = True
                else:
                    assert not stranded  # Don't have both +/- and .
                    strand = "."  # N/A
                    stranded = False

                fields["strand"] = strand

                data[chrom].append(fields)
            except (IndexError, ValueError):
                die("Error parsing fields from feature line:\n\t%s" % line)

    if sort:
        # Sort features by ascending start
        for chrom_features in data.itervalues():
            chrom_features.sort(key=itemgetter("start"))

    return data

def gff2arraydict(filename, gtf=False, sort=None, verbose=True):
    """Parses a gff/gtf file and returns a dict of numpy arrays for each chrom.

    Default behavior is to sort features unless 'gtf' is True

    Returns labels, sequences, where:
      labels: a dict from feature_key -> feature string
      sequences: a dict from seqname -> features (numpy array)
        (features: structured array with start, end, key, and source fields.
         If stranded, a strand field is included.
         If gtf-mode, gene_id and transcript_id fields are included.)
    """
    if sort is None:
        sort = not gtf

    data = defaultdict(list)  # A dictionary-like object
    feature_dict = {}
    source_dict = {}
    last_segment_start = {}  # per seq
    unsorted_seqs = set()
    stranded = None
    with maybe_gzip_open(filename) as infile:
        for datum in read_gff(infile, gtf=gtf):
            try:
                if datum.start < last_segment_start[datum.seqname]:
                    unsorted_seqs.add(datum.seqname)
            except KeyError:
                pass

            feature = str(datum.feature)
            try:  # Lookup feature key
                feature_key = feature_dict[feature]
            except KeyError:  # Map new feature to key
                feature_key = len(feature_dict)
                feature_dict[feature] = feature_key

            source = str(datum.source)
            try:  # Lookup source key
                source_key = source_dict[source]
            except KeyError:  # Map new source to key
                source_key = len(source_dict)
                source_dict[source] = source_key

            try:
                strand = datum.strand
            except AttributeError:
                strand = "."

            if strand == "+" or strand == "-":
                assert stranded is None or stranded
                stranded = True
            else:
                assert not stranded
                stranded = False

            feature = [datum.start, datum.end, feature_key, source_key]
            if stranded:
                feature.append(datum.strand)

            if gtf:
                gene_id = datum.attributes["gene_id"]
                transcript_id = datum.attributes["transcript_id"]
                if len(gene_id) > MAXLEN_ID_DTYPE:
                    raise ValueError("gene_id field was too long (over \
%d characters): %s" % (MAXLEN_ID_DTYPE, gene_id))
                elif len(transcript_id) > MAXLEN_ID_DTYPE:
                    raise ValueError("transcript_id field was too long (over \
%d characters): %s" % (MAXLEN_ID_DTYPE, transcript_id))
                feature.extend([gene_id, transcript_id])

            data[datum.seqname].append(tuple(feature))
            last_segment_start[datum.seqname] = datum.start

    # Create reverse dict for return
    features = inverse_dict(feature_dict)
    sources = inverse_dict(source_dict)

    # convert lists of tuples to array
    dtype = [('start', DTYPE_SEGMENT_START),
             ('end', DTYPE_SEGMENT_END),
             ('key', DTYPE_SEGMENT_KEY),
             ('source_key', DTYPE_SOURCE_KEY)]
    if stranded:
        dtype.append(('strand', DTYPE_STRAND))
    if gtf:
        dtype.extend([('gene_id', DTYPE_GENE_ID),
                      ('transcript_id', DTYPE_TRANSCRIPT_ID)])

    sequences = dict((seq, array(segments, dtype=dtype))
                       for seq, segments in data.iteritems())

    # Sort segments within each chromosome
    for seq, segments in sequences.iteritems():
        if sort and seq in unsorted_seqs:
            if verbose:
                print >>sys.stderr, \
                    "Segments were unsorted relative to %s. Sorting..." % seq,
            segments.sort(order=['start'])
            if verbose:
                print >>sys.stderr, "done"

    return features, sources, sequences

def load_features(filename, verbose=True, *args, **kwargs):
    if verbose:
        print >>sys.stderr, "Loading features from: %s..." % filename

    features, sources, sequences = \
        gff2arraydict(filename, verbose=verbose, *args, **kwargs)

    if verbose:
        print >>sys.stderr, "\tFound %d sources:" % len(sources)
        for source_key, source in sources.iteritems():
            print >>sys.stderr, '\t\t"%s" -> %d' % (source, source_key)

        print >>sys.stderr, "\tFound %d feature classes" % len(features)
        for key, feature in features.iteritems():
            print >>sys.stderr, '\t\t"%s" -> %d' % (feature, key)

        print >>sys.stderr, "done"

    return Features(sequences, features, sources)
