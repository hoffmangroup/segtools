#!/usr/bin/env python
from __future__ import division, with_statement

"""
bed_compare.py

Tools for comparing two bed files.

"""

__version__ = "$Revision: 151 $"

import os
import sys

from .common import die, load_segmentation, SEGMENT_START_COL, SEGMENT_END_COL, SEGMENT_LABEL_KEY_COL

def bases_in_segments(segments):
    """Return the number of bases in a segment array"""
    if segments is None or segments.shape[0] == 0:
        return 0
    else:
        return (segments[:, SEGMENT_END_COL] - \
                    segments[:, SEGMENT_START_COL]).sum()

def edit_distance(bedfile1, bedfile2, quick=False, verbose=False):
    """Given two bed files, prints the edit distance (bp) between them"""
    segmentation1 = load_segmentation(bedfile1, verbose=verbose)
    if segmentation1 is None: die("Failed to load segmentation: %s" % bedfile1)
    segmentation2 = load_segmentation(bedfile2, verbose=verbose)
    if segmentation2 is None: die("Failed to load segmentation: %s" % bedfile2)

    chroms = set(segmentation1.chromosomes.keys() + \
                     segmentation2.chromosomes.keys())
    bases_diff = 0
    bases_same = 0
    bases_missing1 = 0
    bases_missing2 = 0
    for chrom in chroms:
        if verbose:
            print "%s" % chrom

        # If no segments in segmentation 1
        try:
            segs1 = segmentation1.chromosomes[chrom]
        except KeyError:
            segs1 = None


        # If no segments in segmentation 2
        try:
            segs2 = segmentation2.chromosomes[chrom]
        except KeyError:
            segs2 = None

        if segs1 is None and segs2 is None:
            continue
        elif segs1 is None or segs2 is None:
            if segs1 is None:
                bases_missing1 += bases_in_segments(segs2)
            elif segs2 is None:
                bases_missing2 += bases_in_segments(segs1)
            continue

        # Segments in both segmentations
        segs1_iter = iter(segs1)
        start1, end1, label_key1 = segs1_iter.next()
        segs2_iter = iter(segs2)
        start2, end2, label_key2 = segs2_iter.next()
        while True:
            advance1 = False
            advance2 = False
            if start1 < start2:  # move up segment 1
                stop = min(start2, end1)
                bases_missing2 += stop - start1
                start1 = stop
                if end1 <= start2:
                    advance1 = True
            elif start2 < start1:  # move up segment 2
                stop = min(start1, end2)
                bases_missing1 += stop - start2
                start2 = stop
                if end2 <= start1:
                    advance2 = True
            else:  # start1 == start2
                if end1 < end2:
                    bases = end1 - start1
                    stop = end1
                    advance1 = True
                elif end2 < end1:
                    bases = end2 - start2
                    stop = end2
                    advance2 = True
                else:  # Segments match perfectly
                    bases = end2 - start2
                    stop = end2
                    advance1 = True
                    advance2 = True

                if label_key1 == label_key2:
                    bases_same += bases
                else:
                    bases_diff += bases

                # Advance both
                start1 = stop
                start2 = stop

            # Carry out any pending advances
            if advance1:
                try:
                    start1, end1, label_key1 = segs1_iter.next()
                except StopIteration:
                    bases_missing1 += end2 - start2
                    for start2, end2, label_key2 in segs2_iter:
                        bases_missing1 += end2 - start2
                    break

            if advance2:
                try:
                    start2, end2, label_key2 = segs2_iter.next()
                except StopIteration:
                    bases_missing2 += end1 - start1
                    for start1, end1, label_key1 in segs1_iter:
                        bases_missing2 += end1 - start1
                    break

        if quick: break

    # Direct just numbers to stdout
    stderr = ["Bases the same:  ",
              "Bases different: ",
              "Bases missing in %s:\t " % bedfile1,
              "Bases missing in %s:\t " % bedfile2]
    stdout = ["%d" % bases_same,
              "%d" % bases_diff,
              "%d" % bases_missing1,
              "%d" % bases_missing2]

    print >>sys.stderr, "\n===== EDIT DISTANCE ====="
    mesh_output(stderr, stdout)

def mesh_output(messages, values):
    """Print messages to stderr, values tab-delimited to stdout"""
    assert len(messages) == len(values)
    for message, value in zip(messages, values):
        print >>sys.stderr, str(message),
        sys.stderr.flush()
        print "%s\t" % value,
        sys.stdout.flush()
        print >>sys.stderr, ""
        sys.stderr.flush()

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] BEDFILE BEDFILE"
    description = "Compare two BED-formatted segmentations"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-d", "--edit-distance", dest="edit_distance",
                      default=False, action="store_true",
                      help="Measure the base-wise edit distance between"
                      " the two specified segmentations")
    parser.add_option("-q", "--quick", dest="quick",
                      default=False, action="store_true",
                      help="Output results after only one chromosome")
    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print diagnostic messages")

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfiles = args[0:2]

    if options.edit_distance:
        kwargs = {"quick": options.quick,
                  "verbose": options.verbose}
        edit_distance(*bedfiles, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
