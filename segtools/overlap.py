#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
seg_overlap.py:

Evaluates the overlap between two BED files, based upon the spec at:
http://encodewiki.ucsc.edu/EncodeDCC/index.php/Overlap_analysis_tool_specification


Author: Orion Buske <orion.buske@gmail.com>
Date:   August 18, 2009
"""

import math
import os
import sys

from math import ceil
from numpy import invert, logical_or, zeros
from rpy2.robjects import r, numpy2ri

# XXX: Do this without the kludgy constants
from .common import die, get_ordered_labels, image_saver, load_segmentation, make_tabfilename, map_mnemonics, maybe_gzip_open, OutputMasker, r_source, SEGMENT_START_COL, SEGMENT_END_COL, SEGMENT_LABEL_KEY_COL, setup_directory, tab_saver, template_substitute

from .bed import read_native
from .html import find_output_files, save_html_div

# A package-unique, descriptive module name used for filenames, etc
MODULE = "overlap"

NAMEBASE = "%s" % MODULE
FORWARD_NAMEBASE = os.extsep.join([NAMEBASE, "seg"])
BACKWARD_NAMEBASE = os.extsep.join([NAMEBASE, "feature"])
SIGNIFICANCE_NAMEBASE = os.extsep.join([NAMEBASE, "significance"])

HTML_TITLE_BASE = "Overlap statistics"
HTML_TEMPLATE_FILENAME = "overlap_div.tmpl"
SIGNIFICANCE_TEMPLATE_FILENAME = "overlap_significance.tmpl"

ROW_NAME_COL = ""
NONE_COL = "none"
TOTAL_COL = "total"

BY_CHOICES = ["segments", "bases"]
BY_DEFAULT = "segments"
MIDPOINT_CHOICES = ["1", "2", "both"]
SAMPLES_DEFAULT = 1000
REGION_FRACTION_DEFAULT = 0.5
SUBREGION_FRACTION_DEFAULT = 0.5

PNG_SIZE_PER_PANEL = 400  # px
SIGNIFICANCE_PNG_SIZE = 600  # px

def start_R():
    r_source("common.R")
    r_source("overlap.R")

def calc_overlap(subseg, qryseg, quick=False, by=BY_DEFAULT,
                 min_overlap=1, min_overlap_fraction=None):
    # Ensure either min_overlap or _fraction, but not both
    if min_overlap_fraction is None:
        min_overlap = int(min_overlap)
    else:
        min_overlap = None
        min_overlap_fraction = float(min_overlap_fraction)
        assert min_overlap_fraction >= 0 and min_overlap_fraction <= 1

    sub_labels = subseg.labels
    qry_labels = qryseg.labels
    # Make sure reserved column names aren't used
    for col_label in [ROW_NAME_COL, NONE_COL, TOTAL_COL]:
        if col_label in sub_labels.values():
            die("SUBJECTFILE uses reserved group name: %s" % col_label)
        elif col_label in qry_labels.values():
            die("QUERYFILE uses reserved group name: %s" % col_label)

    counts = zeros((len(sub_labels), len(qry_labels)), dtype="int")
    totals = zeros(len(sub_labels), dtype="int")
    nones = zeros(len(sub_labels), dtype="int")
    for chrom in subseg.chromosomes:
        print >>sys.stderr, "\t%s" % chrom
        try:
            qry_segments = qryseg.chromosomes[chrom]
        except KeyError:
            print >>sys.stderr, "\t\tNot in QUERYFILE. Skipping."
            continue

        # Track upper and lower bounds on range of segments that might
        #   be in the range of the current segment (to keep O(n))
        low_qry_segment_iter = iter(qry_segments)
        high_qry_segment_iter = iter(qry_segments)
        low_index = 0
        low_stop = False
        low_qry_segment = low_qry_segment_iter.next()
        high_index = 0
        high_stop = False
        high_qry_segment = high_qry_segment_iter.next()

        # Keep track of the number of segments of each group that are currently
        #   in range of the current segment
        qrygroup_counts = zeros(len(qry_labels), dtype="int")

        for sub_segment in subseg.chromosomes[chrom]:
            substart = sub_segment[SEGMENT_START_COL]
            subend = sub_segment[SEGMENT_END_COL]
            sublabelkey = sub_segment[SEGMENT_LABEL_KEY_COL]
            sublength = subend - substart
            # Compute min-overlap in terms of bases, if conversion required
            if min_overlap is None:
                # Be conservative with ceil
                min_overlap_bp = ceil(min_overlap_fraction * sublength)
            else:
                min_overlap_bp = min_overlap

            # High_index is just above top of range around subsegment
            while not high_stop and \
                    high_qry_segment[SEGMENT_START_COL] <= \
                    sub_segment[SEGMENT_END_COL] - min_overlap:
                try:
                    qrygroup_counts[
                        high_qry_segment[SEGMENT_LABEL_KEY_COL]] += 1
                    high_qry_segment = high_qry_segment_iter.next()
                    high_index += 1
                except StopIteration:
                    high_index = len(qry_segments)
                    high_stop = True
                    break

            # Low_index is bottom of range around subsegment
            while not low_stop and \
                    low_qry_segment[SEGMENT_END_COL] - \
                    min_overlap < sub_segment[SEGMENT_START_COL]:
                try:
                    qrygroup_counts[
                        low_qry_segment[SEGMENT_LABEL_KEY_COL]] -= 1
                    low_qry_segment = low_qry_segment_iter.next()
                    low_index += 1
                except StopIteration:
                    low_index = len(qry_segments)
                    low_stop = True
                    break

            if by == "segments":
                full_score = 1
            elif by == "bases":
                full_score = sublength

            totals[sublabelkey] += full_score

            # Skip processing if there are no overlapping segments
            if high_index - low_index <= 0:
                nones[sublabelkey] += full_score
                continue

            # At least one segment overlaps, so calculate overlap
            range_segments = qry_segments[low_index:high_index, :]
            label_col = range_segments[:, SEGMENT_LABEL_KEY_COL]

            if by == "segments":
                # Add 1 to count for each group currently in range
                groups_in_range = (qrygroup_counts > 0)
                assert groups_in_range.sum() != 0
                counts[sublabelkey] += groups_in_range

            elif by == "bases":
                # Keep track of total covered by all labels
                covered = zeros(sublength, dtype="bool")
                for qrylabelkey in qry_labels:
                    # Look at total covered by single label
                    label_rows = (label_col == qrylabelkey)
                    label_covered = zeros(sublength, dtype="bool")
                    for qry_segment in range_segments[label_rows, :]:
                        qrystart = qry_segment[SEGMENT_START_COL]
                        qryend = qry_segment[SEGMENT_END_COL]
                        # Define bounds of coverage
                        cov_start = min(qrystart - substart, 0)
                        cov_end = max(cov_start + (qryend - qrystart),
                                      sublength)
                        label_covered[cov_start:cov_end]=True

                    # Add the number of bases covered by this segment
                    counts[sublabelkey][qrylabelkey] += label_covered.sum()
                    covered = logical_or(covered, label_covered)

                # See how many bases were never covered by any segment
                nones[sublabelkey] += invert(covered).sum()

        if quick: break
    return (counts, totals, nones)

def split_bed_regions(filename, labels, lengths):
    """Simultaneously splits a bed file into region objects for each label.

    Mostly a direct copy of code in Nathan Boley's GSC
    (encode.py:parse_bed_file)
    Should be ~25 times faster than re-parsing with the 'group' arg for each
    label
    """
    from encode import binary_region, interval, regions

    # Initialize dict of regions
    chrom_dict = {}
    for label in labels.values():
        chrom_dict[label] = regions()
        for chName in lengths.keys():
            length = lengths[chName].length
            chrom_dict[label][chName] = binary_region((), chName, length)

    with maybe_gzip_open(filename) as infile:
        for datum in read_native(infile):
            feature = interval(int(datum.chromStart), int(datum.chromEnd))
            chName = datum.chrom
            name = datum.name
            if not lengths.has_key( chName ):
                continue

            shifted_feature = lengths[chName].intersection(feature) 
            # if there is no intersection, there is nothing else to do
            if shifted_feature == None: continue
            else:
                try:
                    chrom_dict[name][chName].add(shifted_feature)
                # This assumes that, if the region name is not present
                # in the lengths file, then we want to ignore it
                except KeyError:
                    pass

    return chrom_dict

def parse_file_regions(parser, filename, labels, lengths):
    regions = {}
    for key, label in labels.iteritems():
        print >>sys.stderr, "\t%s" % label
        with maybe_gzip_open(filename) as file:
            regions[key] = parser(file, lengths,
                                  group=label,
                                  force_binary=True)
    return regions

def calc_significance(seg_labels, feature_labels, bedfilename,
                      featurefilename, regionfilename,
                      num_samples=SAMPLES_DEFAULT,
                      region_fraction=REGION_FRACTION_DEFAULT,
                      subregion_fraction=SUBREGION_FRACTION_DEFAULT,
                      quick=False, by=BY_DEFAULT):
    """Use Genome Structure Correction tool to measure overlap significance.

    If GSC is not in the user's PYTHONPATH or the user failed to provide a
    region file for GSC, skip this part of validation and move on.

    Returns a nested dictionary mapping label-feature pairs to p_values:
      dict ( seg_label_key -> dict ( feature_label_key -> p_value ) )
    """
    try:
        import block_bootstrap
    except ImportError:
        print >> sys.stderr, "GSC not on PYTHONPATH. \
Skipping significance statistics."
        return None
    else:
        print >>sys.stderr, "Found GSC"

    from block_bootstrap import marginal_bp_overlap_stat, \
        marginal_resample_region_overlap_stat
    block_bootstrap.output = OutputMasker()  # Block GSC stdout output
    from encode import parse_lengths_file

    if regionfilename is None:
        print >> sys.stderr, "A region file must be specified to use GSC. \
Skipping significance statistics."
        return None
    elif not os.path.isfile(regionfilename):
        print >> sys.stderr, "Could not locate region file: %s. \
Skipping significance statistics." % regionfilename
        return None
    else:
        print >>sys.stderr, "Using region file: %s" % regionfilename

    # Parse region file
    print >>sys.stderr, "Loading region file...",
    sys.stdout.flush()
    with maybe_gzip_open(regionfilename) as regionfile:
        lengths = parse_lengths_file(regionfile)
    print >>sys.stderr, "done"

    # Preprocess both input BED files, splitting by labels

    # Parse segmentation BED file, once for each label
    print >>sys.stderr, "Parsing segmentation by label...",
    sys.stdout.flush()
    seg_regions = split_bed_regions(bedfilename, seg_labels, lengths)
    #seg_regions = parse_file_regions(parse_bed_file, bedfilename,
    #                                 seg_labels, lengths)
    print >>sys.stderr, "done"

    # Parse feature BED file, once for each label
    print >>sys.stderr, "Parsing features by label...",
    sys.stdout.flush()
    feature_regions = split_bed_regions(featurefilename,
                                        feature_labels, lengths)
    #feature_regions = parse_file_regions(parse_bed_file, featurefilename,
    #                                     feature_labels, lengths)
    print >>sys.stderr, "done"

    # Calculate p-value for every label-feature pair
    p_values = {}
    print >>sys.stderr, "Calculating p-values:"
    for seg_label_key, seg_label in seg_labels.iteritems():
        print >>sys.stderr, "\t%s" % seg_label
        p_values[seg_label_key] = {}
        for feature_label_key, feature_label in feature_labels.iteritems():
            # Use GSC to calculate p-value
            try:
                sys.stdout = OutputMasker(sys.stdout)  # Kill stdout prints
                if by == "segments":
                    z, p = marginal_resample_region_overlap_stat(
                        seg_regions[seg_label], feature_regions[feature_label],
                        region_fraction, subregion_fraction, num_samples)
                elif by == "bases":
                    z, p = marginal_bp_overlap_stat(
                        seg_regions[seg_label], feature_regions[feature_label],
                        region_fraction, num_samples)
                else:
                    print >> sys.stderr, "Unrecognized mode: %s.\
 Skipping significance testing." % by
                    return None

            except Exception, e:
                print >> sys.stderr, "GSC raised the following exception:\
 %s.\nSkipping significance testing." % e
                return None
            finally:
                sys.stdout = sys.stdout.restore()

            p_values[seg_label_key][feature_label_key] = p

        if quick: break

    return p_values

def make_tab_row(row_labels, row_key, col_labels, col_ordered_keys,
                 data, nones=None, totals=None, formatstr="%d"):
    row = {ROW_NAME_COL : row_labels[row_key]}
    for col_key in col_ordered_keys:
        row[col_labels[col_key]] = formatstr % data[row_key][col_key]

    if nones is not None:
        row[NONE_COL] = nones[row_key]

    if totals is not None:
        row[TOTAL_COL] = totals[row_key]

    return row

## Saves the data to a tab file
def save_tab(dirpath, row_labels, col_labels, counts, totals, nones,
             namebase=NAMEBASE, clobber=False,
             mnemonic_rows=[], mnemonic_cols=[]):
    assert counts is not None and totals is not None and nones is not None

    row_ordered_keys, row_labels = get_ordered_labels(row_labels,
                                                      mnemonic_rows)
    col_ordered_keys, col_labels = get_ordered_labels(col_labels,
                                                      mnemonic_cols)

    # Set up fieldnames based upon QUERYFILE groups
    fieldnames = list(col_labels[col_key] for col_key in col_ordered_keys)
    fieldnames.insert(0, ROW_NAME_COL)
    fieldnames.append(NONE_COL)
    fieldnames.append(TOTAL_COL)
    with tab_saver(dirpath, namebase, fieldnames,
                   clobber=clobber) as count_saver:
        for row_key in row_ordered_keys:
            count_row = make_tab_row(row_labels, row_key, col_labels,
                                     col_ordered_keys, counts,
                                     nones, totals)
            count_saver.writerow(count_row)

## Saves the significance data to a tab file
def save_significance_tab(dirpath, seg_labels, feature_labels, p_values,
             namebase=SIGNIFICANCE_NAMEBASE, clobber=False,
             mnemonics=[]):
    assert p_values is not None

    row_ordered_keys, row_labels = get_ordered_labels(seg_labels, mnemonics)
    col_ordered_keys, col_labels = get_ordered_labels(feature_labels)

    # Set up fieldnames based upon feature groups
    fieldnames = list(col_labels[col_key] for col_key in col_ordered_keys)
    fieldnames.insert(0, ROW_NAME_COL)
    with tab_saver(dirpath, namebase, fieldnames,
                   clobber=clobber) as count_saver:
        for row_key in row_ordered_keys:
            if row_key in p_values:
                count_row = make_tab_row(row_labels, row_key, col_labels,
                                         col_ordered_keys, p_values,
                                         formatstr="%.0e")
                count_saver.writerow(count_row)

def save_plot(dirpath, num_panels, clobber=False, mnemonics=[]):
    start_R()

    tabfilename = make_tabfilename(dirpath, FORWARD_NAMEBASE)
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
                                 dirpath=dirpath,
                                 basename=NAMEBASE,
                                 mnemonics=mnemonics))

def save_significance_plot(dirpath, clobber=False, mnemonics=[]):
    start_R()

    tabfilename = make_tabfilename(dirpath, SIGNIFICANCE_NAMEBASE)
    if not os.path.isfile(tabfilename):
        print >> sys.stderr, "Unable to find tab file: %s. \
Skipping significance plot." % tabfilename
        return

    # Plot data in file
    with image_saver(dirpath, SIGNIFICANCE_NAMEBASE, clobber=clobber,
                     width=SIGNIFICANCE_PNG_SIZE,
                     height=SIGNIFICANCE_PNG_SIZE):
        r.plot(r["plot.overlap.pvalues"](tabfilename, mnemonics=mnemonics))

def save_html(dirpath, bedfilename, featurefilename, by, clobber=False):
    bedfilename = os.path.basename(bedfilename)
    featurebasename = os.path.basename(featurefilename)

    title = "%s (%s)" % (HTML_TITLE_BASE, featurebasename)

    # Insert significance results if they were made
    tabfilename = make_tabfilename(dirpath, SIGNIFICANCE_NAMEBASE)
    if os.path.isfile(tabfilename):
        files = find_output_files(dirpath, SIGNIFICANCE_NAMEBASE)
        significance = template_substitute(
            SIGNIFICANCE_TEMPLATE_FILENAME)(files)
    else:
        significance = ""

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  title=title, forwardtablenamebase=FORWARD_NAMEBASE,
                  backwardtablenamebase=BACKWARD_NAMEBASE,
                  module=MODULE, by=by, significance=significance,
                  bedfilename=bedfilename, featurefilename=featurebasename)

## Package entry point
def validate(bedfilename, featurefilename, dirpath, regionfilename=None,
             clobber=False, quick=False,
             by=BY_DEFAULT, samples=SAMPLES_DEFAULT,
             region_fraction=REGION_FRACTION_DEFAULT,
             subregion_fraction=SUBREGION_FRACTION_DEFAULT,
             min_overlap=1, min_overlap_fraction=None,
             mnemonicfilename=None, replot=False, noplot=False):
    setup_directory(dirpath)
    segmentation = load_segmentation(bedfilename)
    features = load_segmentation(featurefilename)

    assert segmentation is not None
    assert features is not None

    seg_labels = segmentation.labels
    feature_labels = features.labels
    mnemonics = map_mnemonics(seg_labels, mnemonicfilename)

    if not replot:
        # Overlap of segmentation with features (forward)
        print >>sys.stderr, "Measuring overlap of segmentation with features..."
        seg_counts, seg_totals, seg_nones = calc_overlap(
            segmentation, features, by=by, min_overlap=min_overlap,
            min_overlap_fraction=min_overlap_fraction, quick=quick)
        save_tab(dirpath, seg_labels, feature_labels, seg_counts,
                 seg_totals, seg_nones, namebase=FORWARD_NAMEBASE,
                 clobber=clobber, mnemonic_rows=mnemonics)  # Labels on rows

        # Overlap of features with segmentation (backward)
        print >>sys.stderr, "Measuring overlap of features with segmentation..."
        feature_counts, feature_totals, feature_nones=calc_overlap(
            features, segmentation, by=by, min_overlap=min_overlap,
            min_overlap_fraction=min_overlap_fraction, quick=quick)
        save_tab(dirpath, feature_labels, seg_labels, feature_counts,
                 feature_totals, feature_nones, namebase=BACKWARD_NAMEBASE,
                 clobber=clobber, mnemonic_cols=mnemonics)  # Labels on cols

        # GSC significance of forward overlap
        print >>sys.stderr, "Measuring significance of overlap \
of segmentation with features..."
        p_values = calc_significance(seg_labels, feature_labels, bedfilename,
                                     featurefilename, regionfilename,
                                     num_samples=samples,
                                     region_fraction=region_fraction,
                                     subregion_fraction=subregion_fraction,
                                     quick=quick, by=by)
        if p_values is not None:
            save_significance_tab(dirpath, seg_labels, feature_labels,
                                  p_values=p_values, clobber=clobber,
                                  mnemonics=mnemonics)

    if not noplot:
        save_plot(dirpath, num_panels=len(feature_labels),
                  clobber=clobber, mnemonics=mnemonics)
        save_significance_plot(dirpath, clobber=clobber,
                               mnemonics=mnemonics)

    print >>sys.stderr, "Saving HTML div...",
    sys.stdout.flush()  # Necessary to make sure html errors don't clobber print
    save_html(dirpath, bedfilename, featurefilename, by=by, clobber=clobber)
    print >>sys.stderr, "done"

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] BEDFILE FEATUREFILE"
    description = "BEDFILE and FEATUREFILE should both be in BED3+ format \
(gzip'd okay). BEDFILE should correspond to a segmentation. Overlap \
analysis will be performed in both directions (BEDFILE as SUBJECTFILE \
and QUERYFILE). See for full specification: \
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
    group.add_option("--replot", action="store_true",
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate plots")
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
    group.add_option("--mnemonic-file", dest="mnemonicfilename",
                     default=None,
                     help="If specified, labels will be shown using"
                     " mnemonics found in this file")
    group.add_option("-o", "--outdir",
                     dest="outdir", default="%s" % MODULE,
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "GSC Options")
    group.add_option("--region-file", dest="regionfilename",
                     default=None,
                     help="If specified, this file will be used to calculate"
                     " overlap significance using GSC. This must be a BED file")
    group.add_option("-s", "--samples", type="int",
                     dest="samples", default=SAMPLES_DEFAULT,
                     help="The number of samples for GSC to use to estimate"
                     " the significance of the overlap [default: %default]")
    group.add_option("--region-fraction", type="float",
                     dest="region_fraction",
                     default=REGION_FRACTION_DEFAULT,
                     help="The region_fraction tu use with GSC"
                     " [default: %default]")
    group.add_option("--subregion-fraction", type="float",
                     dest="subregion_fraction",
                     default=SUBREGION_FRACTION_DEFAULT,
                     help="The subregion_fraction tu use with GSC"
                     " [default: %default]")
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
    kwargs = {"regionfilename": options.regionfilename,
              "clobber": options.clobber,
              "quick": options.quick,
              "replot": options.replot,
              "noplot": options.noplot}
              "by": options.by,
              "min_overlap": options.min_overlap,
              "samples": options.samples,
              "min_overlap_fraction": options.min_overlap_fraction,
              "region_fraction": options.region_fraction,
              "subregion_fraction": options.subregion_fraction,
              "mnemonicfilename": options.mnemonicfilename,
    validate(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
