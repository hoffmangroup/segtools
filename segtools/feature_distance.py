#!/usr/bin/env python

"""
Given a SEGMENTATION and ANNOTATION file, 1) prints the distance of each
segment to the nearest feature in the ANNOTATION file
(zero if the two overlap) and
2) generates a histogram of these distances.
Distance is the difference between the nearest bases of the segment and the
feature, so if there is one base pair between them, the distance is 2.
Distances can be strand-corrected with respect to stranded features
by specifying --stranded.
"""

from __future__ import division, with_statement


import os
import sys
from time import time

from numpy import NaN, zeros

from . import Annotation, log, Segmentation, RInterface, open_transcript, \
     add_common_options
from .common import FeatureScanner, die, make_tabfilename, \
     setup_directory, tab_saver
from .html import save_html_div
from .version import __version__

MODULE = "feature_distance"
NAMEBASE = MODULE
HTML_TITLE_BASE = "Distance to nearest feature"
HTML_TEMPLATE_FILENAME = "distance_div.tmpl"


PRINT_FIELDS = ["chrom", "start", "end", "name", "distance"]
TAB_FIELDS = ["label", "group", "distance", "count"]
DELIM = "\t"

DEFAULT_BINS = 100
DEFAULT_BIN_WIDTH = 100

R = RInterface(["common.R", "distance.R"])

def print_line(labels, chrom, segment, distance):
    fields = [chrom,
              str(segment['start']),
              str(segment['end']),
              labels[segment['key']],
              str(distance)]
    print DELIM.join(fields)

def save_tab(counts, bins, labels, groups, outdir,
             clobber=False, verbose=True):
    label_keys = sorted(labels.keys())
    group_keys = sorted(groups.keys())
    with tab_saver(outdir, NAMEBASE, TAB_FIELDS, verbose=verbose,
                   clobber=clobber) as saver:
        for label_key in label_keys:
            for group_key in group_keys:
                count_vector = counts[label_key, group_key, :]
                # Only print if there is at least one data point
                if count_vector.sum() == 0: continue
                for bin, count in zip(bins, count_vector):
                    row = dict(zip(TAB_FIELDS,
                                   [labels[label_key], groups[group_key],
                                    str(bin), str(count)]))
                    saver.writerow(row)

def save_plot(outdir, namebase=NAMEBASE, clobber=False, ropts=None,
              verbose=True, mnemonic_file=None, transcriptfile=None):
    R.start(transcriptfile, verbose=verbose)
    # Load data from corresponding tab file
    tabfilename = make_tabfilename(outdir, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.distances", outdir, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber, ropts=ropts)

def save_html(dirpath, feature_file, mnemonic_file=None,
              clobber=False, verbose=True):
    feature_base = os.path.basename(feature_file)

    title = "%s (%s)" % (HTML_TITLE_BASE, feature_base)

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  verbose=verbose, title=title, mnemonic_file=mnemonic_file,
                  module=MODULE, featurefilename=feature_base)

def calc_distance_bin(n_bins, bin_width, distance, stranded):
    if distance == 0:
        return 0

    sign = (distance > 0) * 2 - 1
    bin = int((distance - sign) / bin_width) + sign
    return min(max(bin, -n_bins - 1), n_bins + 1)

def test_calc_distance_bin():
    params = [(10, 5),
              (5, 13),
              (1, 2)]

    for n_bins, bin_width in params:
        stranded = True
        bins = range(n_bins * 2 + 3)
        tests = [(0, 0),
                 (1, 1),
                 (bin_width, 1),
                 (bin_width + 1, 2),
                 (-1, bins[-1]),
                 (-2, bins[-1]),
                 (n_bins * bin_width, n_bins),
                 (n_bins * bin_width + 1, n_bins + 1),
                 (-n_bins * bin_width - 1, n_bins + 2),
                 (-n_bins * bin_width, n_bins + 3)]
        for distance, answer in tests:
            index = calc_distance_bin(n_bins, bin_width, distance, stranded)
            assert bins[index] == answer, """
Failed test:
  bins:      %r
  bin width: %r
  distance:  %r
  found:     %r (%r)
  correct:   %r
""" % (n_bins, bin_width, distance, bins[index], index, answer)

def print_distances(segmentation, annotation,
                    n_bins=DEFAULT_BINS,
                    bin_width=DEFAULT_BIN_WIDTH,
                    verbose=True,
                    stranded=False,
                    print_nearest=False,
                    all_overlapping=False):
    segment_data = segmentation.chromosomes
    feature_data = annotation.chromosomes

    labels = segmentation.labels
    groups = annotation.labels

    print DELIM.join(PRINT_FIELDS)

    bin_distances = range(0, n_bins * bin_width + 1, bin_width) + \
        ["Inf"]
    if stranded:
        bin_distances.extend(["-Inf"] + range(-n_bins * bin_width, 0, bin_width))

    log("Calulating distances...", verbose)
    start_time = time()
    for chrom in segment_data:
        segments = segment_data[chrom]
        try:
            features = feature_data[chrom]
        except KeyError:
            continue

        log("  %s" % chrom, verbose)
        feature_scanner = FeatureScanner(features)

        count_len = n_bins * (stranded + 1) + 2 + stranded
        assert count_len == len(bin_distances)
        counts = zeros((len(labels), len(groups), count_len), dtype="int")
        for segment in segments:
            nearest = feature_scanner.nearest(segment)
            distance = NaN
            if nearest is not None:
                if isinstance(nearest, list):
                    # If we are only printing one, use the first regardless
                    if not all_overlapping:
                        nearest = nearest[0:1]
                else:
                    # Make sure nearest is a list, even if only one element
                    nearest = [nearest]

                # Save multiple distances as a list of them
                for feature in nearest:
                    distance = FeatureScanner.distance(segment, feature)
                    if stranded:
                        try:
                            strand = feature['strand']
                            if strand == ".":
                                strand = None
                        except IndexError:
                            strand = None

                        if strand is None:
                            die("Trying to use strand information but it was"
                                " not found")

                        if (strand == "+" and \
                                segment['end'] < feature['start']) or \
                                (strand == "-" and \
                                     segment['start'] > feature['end']):
                            distance = -distance


                    count_bin = \
                        calc_distance_bin(n_bins, bin_width,
                                          distance, stranded)

                    counts[segment['key'], feature['key'], count_bin] += 1

                    if print_nearest:
                        name = groups[feature['key']]
                        distance_str = "%s %s" % (distance, name)
                    else:
                        distance_str = str(distance)

                    print_line(labels, chrom, segment, distance_str)

    log("Distances calculated in %.1f seconds" % (time() - start_time), verbose)
    return counts, bin_distances

def feature_distance(segment_file, annotation_file,
                     outdir, mnemonic_file=None,
                     n_bins=DEFAULT_BINS, bin_width=DEFAULT_BIN_WIDTH,
                     clobber=False, ropts=None,
                     stranded=False, print_nearest=False,
                     all_overlapping=False,
                     verbose=True, replot=False, noplot=False):
    if not replot:
        setup_directory(outdir)
        segmentation = Segmentation(segment_file, verbose=verbose)
        annotation = Annotation(annotation_file, verbose=verbose)
        labels = segmentation.labels
        groups = annotation.labels
        counts, bins = print_distances(segmentation, annotation,
                                       n_bins=n_bins, bin_width=bin_width,
                                       verbose=verbose,
                                       stranded=stranded,
                                       print_nearest=print_nearest,
                                       all_overlapping=all_overlapping)
        save_tab(counts, bins, labels, groups, outdir,
                 clobber=clobber, verbose=verbose)

    if not noplot:
        with open_transcript(outdir, MODULE) as transcriptfile:
            save_plot(outdir, clobber=clobber, verbose=verbose,
                      mnemonic_file=mnemonic_file, ropts=ropts,
                      transcriptfile=transcriptfile)

    save_html(outdir, annotation_file, mnemonic_file=mnemonic_file,
              clobber=clobber, verbose=verbose)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] SEGMENTATION ANNOTATION"
    version = "%%prog %s" % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=__doc__.strip())

    group = OptionGroup(parser, "Basic options")
    add_common_options(group, ['clobber', 'quiet'])
    group.add_option("-s", "--stranded", dest="stranded",
                     default=False, action="store_true",
                     help="Strand correct features in the ANNOTATION file."
                     " If the feature contains strand information,"
                     " the sign of the distance value is used to portray the"
                     " relationship between the segment and the nearest"
                     " stranded feature. A negative distance value"
                     " indicates that the segment is nearest the 5' end of the"
                     " feature, whereas a positive value indicates the"
                     " segment is nearest the 3' end of the feature.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Print options")
    group.add_option("-p", "--print-nearest", dest="print_nearest",
                     default=False, action="store_true",
                     help="The name of the nearest feature will be printed"
                     " after each distance (with a space separating the"
                     " two) for features from the ANNOTATION file."
                     " If multiple features are equally near (or overlap),"
                     " it is undefined which will be printed")
    group.add_option("-a", "--all-overlapping", dest="all_overlapping",
                     default=False, action="store_true",
                     help="If multiple features in the ANNOTATION file"
                     " overlap a segment, a separate line is printed for each"
                     " overlapping segment-feature pair. This is most"
                     " interesting with --print-nearest. Otherwise,"
                     " the first overlapping segment will be used for"
                     " printing.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Plotting options")
    add_common_options(group, ['noplot', 'replot'])
    group.add_option("-n", "--n-bins", metavar="N", type="int",
                     dest="n_bins", default=DEFAULT_BINS,
                     help="Number of bins to use in histogram for distances"
                     " greater than zero and less than or equal to N*W."
                     " Distances of 0 and greater than N*W are placed in"
                     " additional bins. If --stranded, N more"
                     " bins are included in the negative direction and"
                     " a bin for distances less than -N*W.")
    group.add_option("-w", "--bin-width", metavar="W", type="int",
                     dest="bin_width", default=DEFAULT_BIN_WIDTH,
                     help="Number of bases in each bin."
                     " If --stranded, bins cover distances of"
                     " (-Inf,-N*W), ..., [-W,0), [0], (0,W], ..., (N*W,Inf)."
                     " Otherwise, the bins cover distances of"
                     " [0], (0,W], (W,2W], ..., (N*W,Inf).")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output options")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
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

    feature_distance(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
