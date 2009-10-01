from __future__ import division, with_statement

__version__ = "$Revision$"

"""
common.py:

Assorted utility functions and classes common or useful to most of segtools.


Author: Orion Buske <orion.buske@gmail.com>
"""

import os
import re
import sys

from collections import defaultdict
from contextlib import closing, contextmanager
from functools import partial
from genomedata import Genome
from gzip import open as _gzip_open
from numpy import array, empty, iinfo
from pkg_resources import resource_filename, resource_string
from string import Template
from tabdelim import DictReader, DictWriter
from rpy2.robjects import r, rinterface, numpy2ri
# numpy2ri imported for side-effect

try:
    PKG = __package__  # Python 2.6
except NameError:
    PKG = "segtools"

PKG_R = os.extsep.join([PKG, "R"])
PKG_RESOURCE = os.extsep.join([PKG, "resources"])
    
from .bed import read_native

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

SEGMENT_START_COL = 0
SEGMENT_END_COL = 1
SEGMENT_LABEL_KEY_COL = 2

#FEATURE_FIELDS = ["chrom", "source", "name", "start", "end", "score", "strand"]

r_filename = partial(resource_filename, PKG_R)
template_string = partial(resource_string, PKG_RESOURCE)


## Wrapper for segmentation data object
class Segmentation(object):
    """
    chromosomes: a dict
      key: chromosome name
      val: segments: numpy.ndarray(each row = chromStart, chromEnd, label_key)
           sorted by chromStart
    labels: a dict
      key: int ("label_key")  (a unique id)
      val: printable (str or int)  (what's in the actual BED file)
      
    tracks: a list of track names that were used to obtain the segmentation
    segtool: the tool that was used to obtain this segmentation (e.g. segway)  
    """

    def __init__(self, chromosomes, labels, tracks, segtool):
        self.chromosomes = chromosomes
        self.labels = labels
        self.tracks = tracks
        self.segtool = segtool


class OutputMasker:
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
        
## Wrapper for gff/gtf feature data
# class Feature(object):
#     def __init__(self, line=None, tokens=None):
#         assert (line or tokens) and not (line and tokens)
#         if tokens is None:
#             tokens = line.strip().split("\t")
            
#         self.__dict__ = dict(zip(FEATURE_FIELDS, tokens))
            
#         # Make zero-based, exclusive end:
#         #   http://genome.ucsc.edu/FAQ/FAQformat#format3
#         self.start = int(self.start) - 1
#         self.end = int(self.end)
        
        
## UTILITY FUNCTIONS
        
## Die with error message
def die(msg="Unexpected error"):
    print >> sys.stderr, "\nERROR: %s\n" % msg
    sys.exit(1)
        
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
        die("Output file: %s already exists! Use --clobber to overwrite!" % filename)
        
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

def template_substitute(filename):
    """
    Simplify import resource strings in the package
    """
    return Template(template_string(filename)).substitute
        
@contextmanager
def tab_saver(dirpath, basename, fieldnames, comment=None,
              clobber=False, header=True, verbose=True):
    """
    Saves data to tab file (yeilds DictWriter)
    """
    if verbose:
        print >>sys.stderr, "Saving tab file...",
        
    outfilename = make_tabfilename(dirpath, basename)
    check_clobber(outfilename, clobber)
    with open(outfilename, "w") as outfile:
        if comment is not None and len(comment) > 0:
            if comment != str(comment):  # Not a single string
                for line in comment:
                    print >> outfile, "# %s" % line
            else:
                print >> outfile, "# %s" % comment
        yield DictWriter(outfile, fieldnames,
                         extrasaction="ignore", header=header)
        
    if verbose:
        print >>sys.stderr, "done"

## OLD IMAGE SAVER
## Generator to save an R plot to both png and pdf with only one plot
## Yields to caller to plot, then saves images
# @contextmanager
# def image_saver(dirpath, basename, clobber=False, **kwargs):
#     png_filename = make_filename(dirpath, basename, EXT_PNG)
#     png_slide_filename = make_filename(dirpath,
#                                        "%s%s" % (basename, ".slide"),
#                                        EXT_PNG)
#     pdf_filename = make_filename(dirpath, basename, EXT_PDF)

#     # Make sure there is no silent overwrite
#     filenames = [png_filename, png_slide_filename, pdf_filename]
#     if not clobber:
#         file_exists = list(os.path.isfile(file) for file in filenames)
#         if any(file_exists):
#             die("Output file: %s already exists! Use --clobber to overwrite!" % filenames[file_exists.index(True)])
        
#     # Set up primary device
#     try:
#         r.png(png_filename, **kwargs)
#     except rinterface.RRuntimeError:
#         die('If unable to start PNG device, try `export DISPLAY=":1"` from the shell before re-running validation')
        
#     png_device = r['dev.cur']()
        
#     r['dev.control']("enable")  # Allow copying from device

#     print >>sys.stderr, "\tCreating images..."
#     yield  # Wait for plot

#     try:
#         # Make pdf copy
#         # Dingbats off because otherwise circles displayed as q's
#         r['dev.copy2pdf'](file=pdf_filename, useDingbats=False)
        
#         # Make slide copy
#         r['dev.print.slide'](png_slide_filename)
        
#         # Flush buffer and close plot
#         r['dev.off'](png_device)
        
#         print >>sys.stderr, "\tImages saved!"
#     except:
#         r['dev.off'](png_device)
#         print >> sys.stderr, "ERROR: R crashed. Some images might not be saved."

## CURRENT VERSION OF IMAGE_SAVER
@contextmanager
def image_saver(dirpath, basename, clobber = False, verbose = True,
                downsample = False, **kwargs):
    """
    Generator to save an R plot to both png and pdf with only one plot
    Yields to caller to plot, then saves images
    """
    # Make sure there is no silent overwrite
    filenames = list(make_filename(dirpath, basename, ext)
                     for ext in IMG_EXTS)
    for filename in filenames:
        check_clobber(filename, clobber)

    if verbose:
        print >>sys.stderr, "Creating images...",
        
    try:
        png_filename = make_pngfilename(dirpath, basename)
        r.png(png_filename, **kwargs)
        png_device = r["dev.cur"]()
        r["dev.control"]("enable")
    except rinterface.RRuntimeError:
        die('Image creation failed.\nIf unable to start PNG device, try setting (export/setenv) variable DISPLAY to ":1" from the (bash/c) shell before re-running validation.')
        
    yield  # Wait for plot
    
    # Use R routine to create all the other images (pdf, slide, etc)
    try:
        r["dev.print.images"](basename = basename, dirname = dirpath,
                              downsample = downsample, **kwargs)
        r["dev.off"](png_device)
    except rinterface.RRuntimeError:
        print >> sys.stderr, 'ERROR: Images might not have been saved!'

    if verbose:
        print >>sys.stderr, "done"

def get_bed_metadata(bedfilename):
    with maybe_gzip_open(bedfilename) as ifp:
        return parse_bed_header(ifp.readline())

## Gets track/tool info from BED file header if it exists
def parse_bed_header(line):
    regexp = re.compile('description="(.*) segmentation of (.*)"')

    segtool = "Missing from BED file"
    tracks = ["Missing from BED file"]

    matching = regexp.search(line)
    if matching:
        segtool = matching.group(1)
        tracks = matching.group(2).split(', ')
    return (segtool, tracks)


## Maps label_keys over segments
def map_segments(segments, labels, chrom_size):
    """
    converts a segment mapping into label_keys at every nucleotide position 
        in the chromosome
        
    e.g.  segment labeled 0 at position [0,4) followed by 
          segment labeled 1 at position [4,7) gets converted into:
          0000111
    """  
    segments_dtype = segments.dtype  # MA: the data type of segments
    segments_dtype_max = iinfo(segments_dtype).max # sentinel value		
    # MA: the maximum value supported by the type
    assert segments_dtype_max not in labels.keys()

    res = fill_array(segments_dtype_max, (chrom_size,), segments_dtype)

    # will overwrite in overlapping case
    for segment in segments:
        res[segment[SEGMENT_START_COL]:
                segment[SEGMENT_END_COL]] = segment[SEGMENT_LABEL_KEY_COL]

    return res

## Yields segment and the continuous corresponding to it, for each segment
##   in the chromosome inside of a supercontig
def loop_segments_continuous(chromosome, segmentation, verbose=True):
    try:
        segments = segmentation.chromosomes[chromosome.name]
    except KeyError:
        if verbose:
            print >>sys.stderr, "\t\tskipping: no data"
        raise StopIteration
    
    supercontig_iter = chromosome.itercontinuous()
    supercontig = None
    supercontig_last_start = 0
    segment_i = 0
    num_segments = len(segments)
    for segment in segments:
        segment_i += 1
        seg_start = segment[SEGMENT_START_COL]
        seg_end = segment[SEGMENT_END_COL]  # Exclusive

        while supercontig is None or seg_start >= supercontig.end:
            # Raise StopIteration if out of supercontigs
            supercontig, continuous = supercontig_iter.next()
            # Enforce increasing supercontig indices
            assert supercontig.start >= supercontig_last_start
            supercontig_last_start = supercontig.start

            if verbose:
                print >>sys.stderr, "\n\t\t%s) %d : %d" % (supercontig.name,
                                             supercontig.start,
                                             supercontig.end)

        if verbose and segment_i % 100 == 0:
            print >>sys.stderr, "\r\t\t\tSegment %d / %d" % (segment_i, num_segments),
            sys.stdout.flush()
            
        if seg_end <= supercontig.start:
            continue  # Get next segment

        yield segment, continuous[seg_start:seg_end]

    if verbose:
        print >>sys.stderr, "\r\t\t\tSegment %d / %d" % (segment_i, num_segments)
        

## Returns splice of continuous or sequece from supercontig within segment
## range. Datatype should be "continuous" or "sequence".
## Assumes supercontigs are non-overlapping and sorted by start ascending
## Returns empty array if no supercontig touches the given range
def get_supercontig_splice(chromosome, start, end, datatype="continuous"):
    if datatype not in ("continuous", "sequence"):
        die("Unknown data type in 'get_supercontig_splice!")
        
    # iterate through each supercontig
    saved = array([])
    for supercontig, continuous in chromosome.itercontinuous():
        #print >>sys.stderr, "Examining supercontig [%d, %d]" % (supercontig.start, supercontig.end)

        # Select appropriate data stream
        if datatype == "continuous":
            data = continuous
        else:
            data = supercontig.seq
            
        # XXX Make sure there are no OBOBS in indexing!
        if start >= supercontig.start and end <= supercontig.end:
            # Found supercontig that completely covers range
            return data[start - supercontig.start : end - supercontig.start]
        elif start >= supercontig.start and start < supercontig.end:
            # Overlapping supercontig end
            saved = array(data[start - supercontig.start : len(data) - 1])
        elif end > supercontig.start and end < supercontig.end:
            # Overlapping supercontig start
            data_splice = data[0 : end - supercontig.start]
            # Because supercontigs in sorted order, we can return here
            if saved.shape[0] > 0:
                return concatenate((saved, data_splice))
            else:
                return data_splice
        elif end < supercontig.start:
            # We've passed the appropriate range
            break

    return saved

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
            die("Error: Output directory does not have appropriate permissions!")

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
def map_mnemonics(labels, mnemonicfilename, field="mnemonic"):
    if mnemonicfilename is None:
        return []
    
    label_order, label_data = load_mnemonics(mnemonicfilename)
    str_labels = labels.values()

    mnemonics = []
    # Add mapping for labels in mnemonic file
    for old_label in label_order:
        new_label = label_data[old_label][field]
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
               index    description    mnemonic
      e.g.     4    initiation (strong)    IS
    Returns a tuple of (label_order, label_data)
    
    label_order: an ordered list of string label indices,
      corresponding to the preferred display order in plots
      
    label_data: dict
      key: a string label index
      value: a dict with fields including "mnemonic" and "description",
        where description is a several-word description of the label
        and mnemonic is a few-character identifier
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
            label_index = row["index"]
            label_order.append(label_index)
            label_data[label_index] = row
    
    return (label_order, label_data)

def key_by_field(fieldname):
    def entry_key(entry):
        return entry[fieldname]
    return entry_key

## Parses gff file for features
def load_gff_data(gff_filename):
    '''
    Expects data in GFF format (1-indexed locations):
    CHROM<tab>source<tab>feature<tab>START<tab>END<tab>score<tab>STRAND
    chrom, start, end, and strand are required
    strand can be '+'/'-' or '.', but not both in the same file

    Returns tssdata
    tssdata: a dict
      key: string chromosome name (e.g. "chr13")
      val: a list of dicts (ordered by the start index, ascending)
        key: "start", "end", or "strand"
        val: the associated data item, if it exists
             string for strand
             int (zero-indexed) for start and end (exclusive)
    '''
    data = defaultdict(list)
    stranded = None
    with open(gff_filename) as infile:
        for line in infile:
            tokens = line.strip().split("\t")
            chrom = tokens[0]

            name = tokens[2]
            start = int(tokens[3]) - 1  # Make zero-indexed
            # Make zero-indexed and exclusive (these cancel out)
            end = int(tokens[4])

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
                
            feature = dict(name=name, start=start, end=end, strand=strand)
            data[chrom].append(feature)
            
    # Sort features by ascending start        
    for chrom_features in data.itervalues():
        chrom_features.sort(key=key_by_field("start"))
        
    return data

def load_segmentation(filename, checknames=True):
    """Returns a segmentation object derived from the given BED file

    If the labels in the BED file are all integers, the label_keys will
      be integer representations of those labels (unless checknames is False).
    If all of the labels are non-integers, the label_keys will be
      integers corresponding to the order observed.
    If there are some integer and some string labels, an error will be raised
      (unless checknames is False).
    """
    
    print >>sys.stderr, "Loading segmentation...",
    sys.stdout.flush()
    
    # first get the tracks that were used for this segmentation
    segtool, tracks = get_bed_metadata(filename)

    data = defaultdict(list)  # A dictionary-like object
    label_dict = {}
    # read in as lists of tuples
    with maybe_gzip_open(filename) as infile:
        for datum in read_native(infile):
            label = str(datum.name)
            try:  # Lookup label key
                label_key = label_dict[label]
            except KeyError:  # Map new label to key
                if checknames:
                    try:
                        label_key = int(label)
                    except ValueError:
                        checknames = False
                        if len(label_dict) > 0:  # Previous int labels found
                            die("Found both integer and string labels in BED file: %s" % filename)
                if not checknames:  # Sequential if's important
                    label_key = len(label_dict)
                    
                label_dict[label] = label_key
                
            segment = (datum.chromStart, datum.chromEnd, label_key)
            data[datum.chrom].append(segment)

    # Create reverse dict for return
    labels = dict((val, key) for key, val in label_dict.iteritems())
    print >>sys.stderr, ""
    for key, label in labels.iteritems():
        print >>sys.stderr, "%d -> \"%s\"" % (key, label)
            
    # Sort segments within each chromosome by start index
    for segments in data.itervalues():
        segments.sort()
        
    # convert lists of tuples to array
    chromosomes = dict((chrom, array(segments))
                       for chrom, segments in data.iteritems())

    print >>sys.stderr, "done."
    # wrap in a Segmentation object
    return Segmentation(chromosomes, labels, tracks, segtool)

## Returns a genome object from data in a genomedata directory
def load_genome(genomedatadir):
    if genomedatadir is not None and os.path.isdir(genomedatadir):
        genome = Genome(genomedatadir)
        if not genome or genome is None:
            die("Error: Unable to load genome data from directory: %s" % genomedatadir)
        else:
            return genome
    else:
        die("Error: Could not find genome data directory: %s" % genomedatadir)
