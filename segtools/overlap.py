#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
Evaluates the overlap between two BED files, based upon the spec at:
http://encodewiki.ucsc.edu/EncodeDCC/index.php/Overlap_analysis_tool_specification

Author: Orion Buske <stasis@uw.edu>
Date:   August 18, 2009
"""

import math
import os
import sys

from collections import defaultdict
from math import ceil
from numpy import bincount, cast, iinfo, invert, logical_or, zeros
from rpy2.robjects import r, numpy2ri

from .common import check_clobber, die, get_ordered_labels, image_saver, load_features, load_segmentation, make_tabfilename, map_mnemonics, r_source, setup_directory, SUFFIX_GZ, tab_saver

from .html import save_html_div

# A package-unique, descriptive module name used for filenames, etc
MODULE = "overlap"

NAMEBASE = "%s" % MODULE
HEATMAP_NAMEBASE = os.extsep.join([NAMEBASE, "heatmap"])
SIGNIFICANCE_NAMEBASE = os.extsep.join([NAMEBASE, "significance"])
OVERLAPPING_SEGMENTS_NAMEBASE = os.extsep.join([NAMEBASE, "segments"])
OVERLAPPING_SEGMENTS_FIELDS = ["chrom", "start (zero-indexed)",
                               "end (exclusive)", "group",
                               "[additional fields]"]

HTML_TITLE_BASE = "Overlap statistics"
HTML_TEMPLATE_FILENAME = "overlap_div.tmpl"
SIGNIFICANCE_TEMPLATE_FILENAME = "overlap_significance.tmpl"

NONE_COL = "none"
TOTAL_COL = "total"

BY_CHOICES = ["segments", "bases"]
BY_DEFAULT = "segments"
MIDPOINT_CHOICES = ["1", "2", "both"]
SAMPLES_DEFAULT = 5000
REGION_FRACTION_DEFAULT = 0.2
SUBREGION_FRACTION_DEFAULT = 0.2

PNG_SIZE_PER_PANEL = 400  # px
SIGNIFICANCE_PNG_SIZE = 600  # px
HEATMAP_PNG_SIZE = 600 # px

def start_R():
    r_source("common.R")
    r_source("overlap.R")

def calc_overlap(subseg, qryseg, quick=False, clobber=False, by=BY_DEFAULT,
                 print_segments=False, dirpath=None, verbose=True,
                 min_overlap=1, min_overlap_fraction=None):
    # Ensure either min_overlap or _fraction, but not both
    if min_overlap_fraction is None:
        min_overlap = int(min_overlap)
    else:
        min_overlap = None
        min_overlap_fraction = float(min_overlap_fraction)
        assert min_overlap_fraction >= 0 and min_overlap_fraction <= 1

    if print_segments: assert dirpath is not None

    sub_labels = subseg.labels
    qry_labels = qryseg.labels

    # Set up output files if printing overlapping segments
    if print_segments:
        outfiles = {}
        header = "# %s" % "\t".join(OVERLAPPING_SEGMENTS_FIELDS)
        for sub_label_key, sub_label in sub_labels.iteritems():
            outfilename = os.extsep.join([OVERLAPPING_SEGMENTS_NAMEBASE,
                                          sub_label, "txt"])
            outfilepath = os.path.join(dirpath, outfilename)
            check_clobber(outfilepath, clobber=clobber)
            outfiles[sub_label_key] = open(outfilepath, "w")
            print >>outfiles[sub_label_key], header

    counts = zeros((len(sub_labels), len(qry_labels)), dtype="int")
    totals = zeros(len(sub_labels), dtype="int")
    nones = zeros(len(sub_labels), dtype="int")

    for chrom in subseg.chromosomes:
        if verbose:
            print >>sys.stderr, "\t%s" % chrom

        try:
            qry_segments = qryseg.chromosomes[chrom]
        except KeyError:
            segments = subseg.chromosomes[chrom]
            segment_keys = segments['key']
            # Numpy does not currently support using bincount on unsigned
            # integers, so we need to cast them to signed ints first.
            # To do this safely, we need to make sure the max segment key
            # is below the max signed int value.
            dtype_max = iinfo('int32').max
            assert segment_keys.max() < dtype_max  # Can cast to int32
            segment_keys = cast['int32'](segment_keys)

            # Assumes segment keys are non-negative, consecutive integers
            if by == "segments":
                key_scores = bincount(segment_keys)

            elif by == "bases":
                # Weight each segment by its length
                weights = segments['end'] - segments['start']
                key_scores = bincount(segment_keys, weights=weights)
                key_scores.astype("int")

            totals += key_scores
            nones += key_scores
            continue

        # Track upper and lower bounds on range of segments that might
        #   be in the range of the current segment (to keep O(n))
        qry_segment_iter = iter(qry_segments)
        qry_segment = qry_segment_iter.next()
        qry_segments_in_range = []
        for sub_segment in subseg.chromosomes[chrom]:
            substart = sub_segment['start']
            subend = sub_segment['end']
            sublabelkey = sub_segment['key']
            sublength = subend - substart
            # Compute min-overlap in terms of bases, if conversion required
            if min_overlap is None:
                # Be conservative with ceil
                min_overlap_bp = ceil(min_overlap_fraction * sublength)
            else:
                min_overlap_bp = min_overlap

            if by == "segments":
                full_score = 1
            elif by == "bases":
                full_score = sublength

            # Add subject segment to the total count
            totals[sublabelkey] += full_score

            # Remove from list any qry_segments that are now too low
            i = 0
            while i < len(qry_segments_in_range):
                segment = qry_segments_in_range[i]
                if segment['end'] - min_overlap_bp < substart:
                    del qry_segments_in_range[i]
                else:
                    i += 1

            # Advance qry_segment pointer to past sub_segment, updating list
            while qry_segment is not None and \
                    qry_segment['start'] <= subend - min_overlap_bp:
                if qry_segment['end'] - min_overlap_bp >= substart:
                    # qry_segment overlaps with sub_segment
                    qry_segments_in_range.append(qry_segment)
                try:
                    qry_segment = qry_segment_iter.next()
                except StopIteration:
                    qry_segment = None

            # Skip processing if there aren't any segments in range
            if len(qry_segments_in_range) == 0:
                nones[sublabelkey] += full_score
                continue

            # Scan list for subset that actually overlap current segment
            overlapping_segments = []
            for segment in qry_segments_in_range:
                if segment['start'] <= subend - min_overlap_bp:
                    assert segment['end'] - min_overlap_bp >= substart
                    overlapping_segments.append(segment)

            # Skip processing if there are no overlapping segments
            if len(overlapping_segments) == 0:
                nones[sublabelkey] += full_score
                continue

            if print_segments:
                for segment in overlapping_segments:
                    values = [chrom,
                              segment['start'],
                              segment['end'],
                              qry_labels[segment['key']]]
                    # Add a source if there is one
                    try:
                        values.append(qryseg.sources[segment['source_key']])
                    except AttributeError, IndexError:
                        pass
                    # Add any other data in the segment
                    try:
                        values.extend(tuple(segment)[4:])
                    except IndexError:
                        pass
                    values = [str(val) for val in values]
                    print >>outfiles[sublabelkey], "\t".join(values)

            # Organize overlapping_segments by qrylabelkey
            label_overlaps = defaultdict(list)  # Per qrylabelkey
            for segment in overlapping_segments:
                label_overlaps[segment['key']].append(segment)
            label_overlaps = dict(label_overlaps)  # Remove defaultdict

            if by == "segments":
                # Add 1 to count for each group that overlaps at least
                # one segment
                for qrylabelkey in label_overlaps:
                    counts[sublabelkey, qrylabelkey] += 1
            elif by == "bases":
                # Keep track of total covered by any labels
                covered = zeros(sublength, dtype="bool")
                for qrylabelkey, segments in label_overlaps.iteritems():
                    # Look at total covered by single label
                    label_covered = zeros(sublength, dtype="bool")
                    for segment in segments:
                        qrystart = segment['start']
                        qryend = segment['end']
                        qrylabelkey = segment['key']
                        # Define bounds of coverage
                        cov_start = max(qrystart - substart, 0)
                        cov_end = min(qryend - substart, sublength)
                        label_covered[cov_start:cov_end] = True

                    # Add the number of bases covered by this segment
                    counts[sublabelkey, qrylabelkey] += label_covered.sum()
                    covered = logical_or(covered, label_covered)

                # See how many bases were never covered by any segment
                nones[sublabelkey] += invert(covered).sum()

        if quick: break

    if print_segments:
        for outfile in outfiles.itervalues():
            outfile.close()

    return (counts, totals, nones)

def make_tab_row(col_indices, data, none, total):
    row = [data[col_i] for col_i in col_indices]
    row.extend([none, total])
    return row

## Saves the data to a tab file
def save_tab(dirpath, row_labels, col_labels, counts, totals, nones,
             namebase=NAMEBASE, clobber=False):
    assert counts is not None and totals is not None and nones is not None

    row_label_keys, row_labels = get_ordered_labels(row_labels)
    col_label_keys, col_labels = get_ordered_labels(col_labels)
    colnames = [col_labels[label_key] for label_key in col_label_keys]
    with tab_saver(dirpath, namebase, clobber=clobber) as count_saver:
        header = [""] + colnames + [NONE_COL, TOTAL_COL]
        count_saver.writerow(header)
        for row_label_key in row_label_keys:
            row = make_tab_row(col_label_keys, counts[row_label_key],
                               nones[row_label_key], totals[row_label_key])
            row.insert(0, row_labels[row_label_key])
            count_saver.writerow(row)

def save_plot(dirpath, num_panels, clobber=False,
              row_mnemonics=[], col_mnemonics=[]):
    start_R()

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    panels_sqrt = math.sqrt(num_panels)
    width = math.ceil(panels_sqrt)
    height = math.floor(panels_sqrt)

    # Plot data in file
    with image_saver(dirpath, NAMEBASE, clobber=clobber,
                     width=PNG_SIZE_PER_PANEL * width,
                     height=PNG_SIZE_PER_PANEL * height):
        r.plot(r["plot.overlap"](tabfilename,
                                 mnemonics=row_mnemonics,
                                 col_mnemonics=col_mnemonics))

def save_heatmap_plot(dirpath, clobber=False,
                      row_mnemonics=[], col_mnemonics=[]):
    start_R()

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    # Plot data in file
    with image_saver(dirpath, HEATMAP_NAMEBASE, clobber=clobber,
                     width=HEATMAP_PNG_SIZE,
                     height=HEATMAP_PNG_SIZE):
        r.plot(r["plot.overlap.heatmap"](tabfilename,
                                         mnemonics=row_mnemonics,
                                         col_mnemonics=col_mnemonics))

def save_html(dirpath, bedfilename, featurefilename, by, clobber=False):
    bedfilename = os.path.basename(bedfilename)
    featurebasename = os.path.basename(featurefilename)
    extra_namebases = {"heatmap": HEATMAP_NAMEBASE}

    title = "%s (%s)" % (HTML_TITLE_BASE, featurebasename)

    significance = ""
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  title=title, tablenamebase=NAMEBASE,
                  extra_namebases = extra_namebases,
                  module=MODULE, by=by, significance=significance,
                  bedfilename=bedfilename, featurefilename=featurebasename)

def is_file_type(filename, ext):
    """Return True if the filename is of the given extension type (e.g. 'txt')

    Allows g-zipping

    """
    base = os.path.basename(filename)
    if base.endswith(SUFFIX_GZ):
        base = base[:-3]
    return base.endswith("." + ext)

## Package entry point
def validate(bedfilename, featurefilename, dirpath, regionfilename=None,
             clobber=False, quick=False, print_segments=False,
             by=BY_DEFAULT, samples=SAMPLES_DEFAULT,
             region_fraction=REGION_FRACTION_DEFAULT,
             subregion_fraction=SUBREGION_FRACTION_DEFAULT,
             min_overlap=1, min_overlap_fraction=None,
             mnemonic_filename=None, feature_mnemonic_filename=None,
             replot=False, noplot=False, verbose=True):
    setup_directory(dirpath)
    segmentation = load_segmentation(bedfilename)

    if is_file_type(featurefilename, 'bed'):
        features = load_segmentation(featurefilename, verbose=verbose)
    elif is_file_type(featurefilename, 'gff'):
        features = load_features(featurefilename, verbose=verbose)
    elif is_file_type(featurefilename, 'gtf'):
        features = load_features(featurefilename, gtf=True,
                                 sort=True, verbose=verbose)
    else:
        raise NotImplementedError("Only bed, gff, and gtf files are supported \
for FEATUREFILE. If the file is in one of these formats, please use the \
appropriate extension")

    assert segmentation is not None
    assert features is not None

    seg_labels = segmentation.labels
    feature_labels = features.labels
    mnemonics = map_mnemonics(seg_labels, mnemonic_filename)
    feature_mnemonics = map_mnemonics(feature_labels,
                                      feature_mnemonic_filename)
    if not replot:
        # Overlap of segmentation with features
        print >>sys.stderr, "Measuring overlap..."
        counts, nones, totals = \
            calc_overlap(segmentation, features, clobber=clobber,
                         by=by, min_overlap=min_overlap,
                         min_overlap_fraction=min_overlap_fraction,
                         print_segments=print_segments,
                         quick=quick, dirpath=dirpath, verbose=verbose)

        save_tab(dirpath, seg_labels, feature_labels,
                 counts, nones, totals, clobber=clobber)

    if not noplot:
        save_plot(dirpath, num_panels=len(feature_labels),
                  clobber=clobber, row_mnemonics=mnemonics,
                  col_mnemonics=feature_mnemonics)
        save_heatmap_plot(dirpath, clobber=clobber,
                          row_mnemonics=mnemonics,
                          col_mnemonics=feature_mnemonics)

    print >>sys.stderr, "Saving HTML div...",
    sys.stdout.flush()  # Necessary to make sure html errors don't clobber print
    save_html(dirpath, bedfilename, featurefilename, by=by, clobber=clobber)
    print >>sys.stderr, "done"

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] BEDFILE FEATUREFILE"
    description = "BEDFILE must be in BED4+ format (name column used as \
grouping variable). FEATUREFILE should be in BED4+ format or GFF format \
(feature column used as grouping variable). Results summarize the overlap \
of BEDFILE groups with FEATUREFILE groups. The symmetric analysis can \
be performed (if both input files are in BED4+ format) \
by rerunning the program with the input file arguments swapped \
(and a different output directory). A rough specification can be found here: \
http://encodewiki.ucsc.edu/EncodeDCC/index.php/\
Overlap_analysis_tool_specification"

    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Flags")
    group.add_option("--clobber", action="store_true",
                     dest="clobber", default=False,
                     help="Overwrite existing output files if the specified"
                     " directory already exists.")
    group.add_option("--quick", action="store_true",
                     dest="quick", default=False,
                     help="Compute values only for one chromosome.")
    group.add_option("-v", "--verbose", action="store_true",
                     dest="verbose", default=False,
                     help="Print status and diagnostic messages.")
    group.add_option("--replot", action="store_true",
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate plots")
    group.add_option("-p", "--print-segments", action="store_true",
                     dest="print_segments", default=False,
                     help=("For each group"
                     " in the BEDFILE, a separate output file will be"
                     " created that contains a list of all the segments that"
                     " the group was found to overlap with. Output files"
                     " are named %s.X.txt, where X is the name"
                     " of the BEDFILE group.") % OVERLAPPING_SEGMENTS_NAMEBASE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Parameters")
    group.add_option("-b", "--by", choices=BY_CHOICES,
                     dest="by", type="choice", default=BY_DEFAULT,
                     help="One of: "+str(BY_CHOICES)+", which determines the"
                     " definition of overlap. @segments: The value"
                     " associated with two features overlapping will be 1 if"
                     " they overlap, and 0 otherwise. @bases: The value"
                     " associated with two features overlapping will be"
                     " number of base pairs which they overlap."
                     " [default: %default]")
    group.add_option("--midpoint-only", choices=MIDPOINT_CHOICES,
                     dest="midpoint", type="choice", default=None,
                     help="For the specified file (1, 2, or both), use only"
                     "the midpoint of each feature instead of the entire"
                     " width.")
    group.add_option("-m", "--min-overlap", type="int",
                     dest="min_overlap", default=1,
                     help="The minimum number of base pairs that two"
                     " features must overlap for them to be classified as"
                     " overlapping. This integer can be either positive"
                     " (features overlap only if they share at least this"
                     " many bases) or negative (features overlap if there"
                     " are no more than this many bases between them). Both"
                     " a negative min-overlap and --by=bases cannot be"
                     " specified together. [default: %default]")
    group.add_option("--min-overlap-fraction", type="float",
                     dest="min_overlap_fraction", default=None,
                     help="The minimum fraction of the base pairs in the"
                     " subject feature that overlap with the query feature"
                     " in order to be counted as overlapping. Overrides"
                     "--min-overlap.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Files")
    group.add_option("--mnemonic-file", dest="mnemonic_filename",
                     default=None,
                     help="If specified, BEDFILE groups will be shown using"
                     " mnemonics found in this file")
    group.add_option("--feature-mnemonic-file",
                     dest="feature_mnemonic_filename", default=None,
                     help="If specified, FEATUREFILE groups will be shown"
                     " using mnemonics found in this file")
    group.add_option("-o", "--outdir",
                     dest="outdir", default="%s" % MODULE,
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    mof = options.min_overlap_fraction
    if mof is not None and (mof < 0 or mof > 1):
        parser.error("Min-overlap-fraction: %.3f out of range: [0, 1]" % mof)

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    args = [args[0], args[1], options.outdir]
    kwargs = {"clobber": options.clobber,
              "quick": options.quick,
              "replot": options.replot,
              "noplot": options.noplot,
              "print_segments": options.print_segments,
              "by": options.by,
              "min_overlap": options.min_overlap,
              "min_overlap_fraction": options.min_overlap_fraction,
              "mnemonic_filename": options.mnemonic_filename,
              "feature_mnemonic_filename": options.feature_mnemonic_filename,
              "verbose": options.verbose}
    validate(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
