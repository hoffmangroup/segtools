#!/usr/bin/env python

"""
Compare two segmentations
"""

from __future__ import division, with_statement

import sys

from . import log, Segmentation, add_common_options
from .version import __version__

def bases_in_segments(segments):
    """Return the number of bases in a segment array"""
    if segments is None or len(segments) == 0:
        return 0
    else:
        return (segments['end'] - segments['start']).sum()

def edit_distance(bedfile1, bedfile2, quick=False, verbose=True):
    """Given two segmentations, returns the edit distance (bp) between them"""
    segmentations = [Segmentation(bedfile, verbose=verbose)
                     for bedfile in [bedfile1, bedfile2]]

    for segmentation in segmentations:
        for segments in segmentation.chromosomes.itervalues():
            if not (segments['start'][1:] >= segments['end'][:-1]).all():
                raise SyntaxError("Overlapping segments found in: %s" % \
                                      segmentation.name)

    chroms = set()
    for segmentation in segmentations:
        chroms.update(segmentation.chromosomes.keys())

    labels = [segmentation.labels for segmentation in segmentations]

    bases_diff = 0
    bases_same = 0
    bases_missing1 = 0
    bases_missing2 = 0
    for chrom in chroms:
        log("%s" % chrom, verbose)

        # If no segments in segmentation 1
        try:
            segs1 = segmentations[0].chromosomes[chrom]
        except KeyError:
            segs1 = None


        # If no segments in segmentation 2
        try:
            segs2 = segmentations[1].chromosomes[chrom]
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
        start1, end1, key1 = segs1_iter.next()
        segs2_iter = iter(segs2)
        start2, end2, key2 = segs2_iter.next()
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

                if labels[0][key1] == labels[1][key2]:
                    bases_same += bases
                else:
                    bases_diff += bases

                # Advance both
                start1 = stop
                start2 = stop

            # Carry out any pending advances
            if advance1:
                try:
                    start1, end1, key1 = segs1_iter.next()
                except StopIteration:
                    bases_missing1 += end2 - start2
                    for start2, end2, key2 in segs2_iter:
                        bases_missing1 += end2 - start2
                    break

            if advance2:
                try:
                    start2, end2, key2 = segs2_iter.next()
                except StopIteration:
                    bases_missing2 += end1 - start1
                    for start1, end1, key1 in segs1_iter:
                        bases_missing2 += end1 - start1
                    break

        if quick: break

    return bases_same, bases_diff, bases_missing1, bases_missing2

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTATION SEGMENTATION"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=__doc__.strip())

    parser.add_option("-d", "--edit-distance", dest="edit_distance",
                      default=False, action="store_true",
                      help="Measure the base-wise edit distance between"
                      " the two specified segmentations")
    add_common_options(parser, ['quick', 'quiet'])
    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    return (options, args)

def print_edit_distance(*args):
    labels = ["bases the same",
              "bases different",
              "bases missing in first file",
              "bases missing in second file"]
    values = [int(val) for val in args]
    log("\n===== EDIT DISTANCE =====")
    for label, value in zip(labels, values):
        log("%s: \t%s" % (label, value))


## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfiles = args[0:2]
    if options.edit_distance:
        kwargs = dict(options.__dict__)
        results = edit_distance(*bedfiles, **kwargs)
        print_edit_distance(*results)
    else:
        log("No actions were specified")

if __name__ == "__main__":
    sys.exit(main())
