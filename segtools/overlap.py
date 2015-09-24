#!/usr/bin/env python
from __future__ import division, with_statement

"""
Evaluates the overlap between two BED files, based upon the spec at:
http://encodewiki.ucsc.edu/EncodeDCC/index.php/Overlap_analysis_tool_specification
"""

# Author: Orion Buske <stasis@uw.edu>
# Date:   August 18, 2009


import os
import sys

from collections import defaultdict
from numpy import bincount, cast, iinfo, invert, logical_or, zeros

from . import Annotation, log, open_transcript, Segmentation, RInterface, \
     add_common_options
from .common import check_clobber, die, get_ordered_labels, make_tabfilename, \
     setup_directory, SUFFIX_GZ, tab_saver
from .html import save_html_div
from .version import __version__

# A package-unique, descriptive module name used for filenames, etc
MODULE = "overlap"

NAMEBASE = "%s" % MODULE
PERFORMANCE_NAMEBASE = os.extsep.join([NAMEBASE, "performance"])

OVERLAPPING_SEGMENTS_NAMEBASE = os.extsep.join([NAMEBASE, "segments"])
OVERLAPPING_SEGMENTS_FIELDS = ["chrom", "start (zero-indexed)",
                               "end (exclusive)", "group",
                               "[additional fields]"]

HTML_TITLE_BASE = "Overlap statistics"
HTML_TEMPLATE_FILENAME = "overlap_div.tmpl"
SIGNIFICANCE_TEMPLATE_FILENAME = "overlap_significance.tmpl"

NONE_COL = "none"
TOTAL_COL = "total"

MODE_CHOICES = ["segments", "bases"]
MODE_DEFAULT = "bases"
MIDPOINT_CHOICES = ["1", "2", "both"]
SAMPLES_DEFAULT = 5000
REGION_FRACTION_DEFAULT = 0.2
SUBREGION_FRACTION_DEFAULT = 0.2

PNG_SIZE_PER_PANEL = 400  # px
SIGNIFICANCE_PNG_SIZE = 600  # px
HEATMAP_PNG_SIZE = 600 # px

R = RInterface(["common.R", "overlap.R"])

def calc_overlap(subseg, qryseg, quick=False, clobber=False, mode=MODE_DEFAULT,
                 print_segments=False, dirpath=None, verbose=True,
                 min_overlap=1):
    min_overlap = int(min_overlap)

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
            outfiles[sub_label_key].write("%s\n" % header)

    counts = zeros((len(sub_labels), len(qry_labels)), dtype="int64")
    totals = zeros(len(sub_labels), dtype="int64")
    nones = zeros(len(sub_labels), dtype="int64")

    for chrom in subseg.chromosomes:
        log("\t%s" % chrom, verbose)

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
            if mode == "segments":
                key_scores = bincount(segment_keys)

            elif mode == "bases":
                # Weight each segment by its length
                weights = segments['end'] - segments['start']
                key_scores = bincount(segment_keys, weights=weights)
                key_scores.astype("int")
            else:
                raise NotImplementedError("Unknown mode: %r" % mode)

            num_scores = key_scores.shape[0]
            totals[0:num_scores] += key_scores
            nones[0:num_scores] += key_scores
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
            min_overlap_bp = min_overlap

            if mode == "segments":
                full_score = 1
            elif mode == "bases":
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
                    except (AttributeError, IndexError):
                        pass
                    # Add any other data in the segment
                    try:
                        values.extend(tuple(segment)[4:])
                    except IndexError:
                        pass
                    values = [str(val) for val in values]
                    outfiles[sublabelkey].write("%s\n" % "\t".join(values))

            # Organize overlapping_segments by qrylabelkey
            label_overlaps = defaultdict(list)  # Per qrylabelkey
            for segment in overlapping_segments:
                label_overlaps[segment['key']].append(segment)
            label_overlaps = dict(label_overlaps)  # Remove defaultdict

            if mode == "segments":
                # Add 1 to count for each group that overlaps at least
                # one segment
                for qrylabelkey in label_overlaps:
                    counts[sublabelkey, qrylabelkey] += 1
            elif mode == "bases":
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
def save_tab(dirpath, row_labels, col_labels, counts, totals, nones, mode,
             namebase=NAMEBASE, clobber=False, verbose=True):
    assert counts is not None and totals is not None and nones is not None

    row_label_keys, row_labels = get_ordered_labels(row_labels)
    col_label_keys, col_labels = get_ordered_labels(col_labels)
    colnames = [col_labels[label_key] for label_key in col_label_keys]
    metadata = {"mode": mode}
    with tab_saver(dirpath, namebase, clobber=clobber, verbose=verbose,
                   metadata=metadata) as count_saver:
        header = [""] + colnames + [NONE_COL, TOTAL_COL]
        count_saver.writerow(header)
        for row_label_key in row_label_keys:
            row = make_tab_row(col_label_keys, counts[row_label_key],
                               nones[row_label_key], totals[row_label_key])
            row.insert(0, row_labels[row_label_key])
            count_saver.writerow(row)

def save_plot(dirpath, namebase=NAMEBASE, clobber=False, verbose=True,
              row_mnemonic_file=None, col_mnemonic_file=None, ropts=None,
              cluster=False, max_contrast=False, transcriptfile=None):
    R.start(transcriptfile, verbose=verbose)

    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.overlap.heatmap", dirpath, namebase, tabfilename,
           mnemonic_file=row_mnemonic_file,
           col_mnemonic_file=col_mnemonic_file,
           clobber=clobber, cluster=cluster,
           max_contrast=max_contrast, ropts=ropts)

def save_performance_plot(dirpath, namebase=PERFORMANCE_NAMEBASE,
                          clobber=False, verbose=True,
                          row_mnemonic_file=None, col_mnemonic_file=None,
                          transcriptfile=None, ropts=None):
    R.start(transcriptfile, verbose=verbose)

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.overlap.performance", dirpath, namebase, tabfilename,
           mnemonic_file=row_mnemonic_file,
           col_mnemonic_file=col_mnemonic_file,
           clobber=clobber, ropts=ropts)


def save_html(dirpath, bedfilename, featurefilename, mode,
              mnemonicfile=None, clobber=False, verbose=True):
    bedfilename = os.path.basename(bedfilename)
    featurebasename = os.path.basename(featurefilename)
    extra_namebases = {"performance": PERFORMANCE_NAMEBASE}

    title = "%s (%s)" % (HTML_TITLE_BASE, featurebasename)

    significance = ""
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  title=title, tables={"":(PERFORMANCE_NAMEBASE, "exact")},
                  mnemonicfile=mnemonicfile, verbose=verbose,
                  extra_namebases=extra_namebases,
                  module=MODULE, by=mode, significance=significance,
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
def overlap(bedfilename, featurefilename, dirpath, regionfilename=None,
            clobber=False, quick=False, print_segments=False,
            mode=MODE_DEFAULT, samples=SAMPLES_DEFAULT,
            region_fraction=REGION_FRACTION_DEFAULT,
            subregion_fraction=SUBREGION_FRACTION_DEFAULT, min_overlap=1,
            mnemonic_file=None, feature_mnemonic_file=None,
            replot=False, noplot=False, cluster=False, max_contrast=False,
            verbose=True, ropts=None):
    if not replot:
        setup_directory(dirpath)

        segmentation = Segmentation(bedfilename, verbose=verbose)
        features = Annotation(featurefilename, verbose=verbose)

        seg_labels = segmentation.labels
        feature_labels = features.labels

        # Overlap of segmentation with features
        log("Measuring overlap:", verbose)

        counts, nones, totals = \
            calc_overlap(segmentation, features, clobber=clobber,
                         mode=mode, min_overlap=min_overlap,
                         print_segments=print_segments,
                         quick=quick, dirpath=dirpath, verbose=verbose)

        save_tab(dirpath, seg_labels, feature_labels,
                 counts, nones, totals, mode=mode,
                 clobber=clobber, verbose=verbose)

    if not noplot:
        with open_transcript(dirpath, MODULE, True) as transcriptfile:
            # Performance plot is only valid for bases
            if mode == "bases":
                save_performance_plot(\
                    dirpath, clobber=clobber, verbose=verbose,
                    row_mnemonic_file=mnemonic_file, ropts=ropts,
                    col_mnemonic_file=feature_mnemonic_file,
                    transcriptfile=transcriptfile)
            else:
                log("Not in base-mode, so skipping performance plot", verbose)

            save_plot(dirpath, clobber=clobber, cluster=cluster,
                      verbose=verbose, row_mnemonic_file=mnemonic_file,
                      col_mnemonic_file=feature_mnemonic_file, ropts=ropts,
                      max_contrast=max_contrast, transcriptfile=transcriptfile)

    save_html(dirpath, bedfilename, featurefilename, mode=mode,
              mnemonicfile=mnemonic_file, clobber=clobber, verbose=verbose)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] SEGMENTATION ANNOTATION"
    description = "SEGMENTATION and ANNOTATION files should be in BED, GFF, \
or GTF format (grouped on 'name'/'feature' columns). \
Results summarize the overlap \
of SEGMENTATION groups with ANNOTATION groups. The symmetric analysis can \
be performed by rerunning the program with the input file arguments swapped \
(and a different output directory)."

    # A rough specification can be found here: http://encodewiki.ucsc.edu/EncodeDCC/index.php/Overlap_analysis_tool_specification

    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'quick', 'replot',
                               'noplot'])
    group.add_option("--cluster", action="store_true",
                     dest="cluster", default=False,
                     help="Cluster rows and columns in heat map plot")
    group.add_option("-p", "--print-segments", action="store_true",
                     dest="print_segments", default=False,
                     help=("For each group"
                     " in the SEGMENTATION, a separate output file will be"
                     " created that contains a list of all the segments that"
                     " the group was found to overlap with. Output files"
                     " are named %s.X.txt, where X is the name"
                     " of the SEGMENTATION group.")
                     % OVERLAPPING_SEGMENTS_NAMEBASE)
    group.add_option("--max-contrast", action="store_true",
                     dest="max_contrast", default=False,
                     help="Saturate color range instead of having it go from"
                     " 0 to 1")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Parameters")
    group.add_option("-b", "--by", choices=MODE_CHOICES, metavar="MODE",
                     dest="mode", type="choice", default=MODE_DEFAULT,
                     help="One of: "+str(MODE_CHOICES)+", which determines the"
                     " definition of overlap. @segments: The value"
                     " associated with two features overlapping will be 1 if"
                     " they overlap, and 0 otherwise. @bases: The value"
                     " associated with two features overlapping will be"
                     " number of base pairs which they overlap."
                     " [default: %default]")
    group.add_option("--min-overlap", type="int",
                     dest="min_overlap", default=1, metavar="N",
                     help="The minimum number of base pairs that two"
                     " features must overlap for them to be classified as"
                     " overlapping. This integer can be either positive"
                     " (features overlap only if they share at least this"
                     " many bases) or negative (features overlap if there"
                     " are no more than this many bases between them). Both"
                     " a negative min-overlap and --by=bases cannot be"
                     " specified together. [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Files")
    add_common_options(group, ['mnemonic_file'])
    group.add_option("--feature-mnemonic-file", metavar="FILE",
                     dest="feature_mnemonic_file", default=None,
                     help="If specified, ANNOTATION groups will be shown"
                     " using mnemonics found in FILE.")
    add_common_options(group, ['outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    kwargs = dict(options.__dict__)
    args.append(kwargs.pop('outdir'))
    overlap(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
