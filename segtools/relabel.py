#!/bin/env python

"""
Generates a new segmentation by first relabeling the segments in a given
SEGMENTATION according to the label mappings in MNEMONICFILE and then
merging overlapping segments with the same final label.
Outputs the new segmentation in bed format to stdout (-o to change).
"""

# Author: Orion Buske <stasis@uw.edu>
# Date:   05.16.2011

from __future__ import division, with_statement

__version__ = "$Revision$"

import os
import sys

from . import log, Segmentation, die, add_common_options
from .common import map_mnemonics, get_ordered_labels


def print_bed(segments, filename=None):
    """Print segments in BED format to file (or stdout if None)"""
    if filename is None:
        outfile = sys.stdout
    else:
        if os.path.isfile(filename):
            log("Warning: overwriting file: %s" % filename)

        outfile = open(filename, "w")

    for segment in segments:
        outfile.write("%s\n" % "\t".join([str(val) for val in segment]))

    if outfile is not sys.stdout:
        outfile.close()

def relabel(segfile, mnemonicfile, outfile=None, verbose=True):
    assert os.path.isfile(segfile)
    segmentation = Segmentation(segfile, verbose=verbose)

    labels = segmentation.labels
    mnemonics = map_mnemonics(labels, mnemonicfile)
    ordered_labels, mnemonics = get_ordered_labels(labels, mnemonics)

    print >>sys.stderr, "Found labels:\n%s" % labels

    print >>sys.stderr, "Found mnemonics:\n%s" % mnemonics

    if outfile is None:
        out = sys.stdout
    else:
        out = open(outfile, "w")


    for chrom, segments in segmentation.chromosomes.iteritems():
        def print_segment(segment, label):
            tokens = [chrom, segment['start'], segment['end'],
                      mnemonics[segment['key']]]
            try:
                tokens.extend([0, segment['strand']])
            except IndexError:
                pass
            out.write('\t'.join([str(token) for token in tokens]))
            out.write('\n')

        prev_segment = segments[0]
        prev_label = mnemonics[prev_segment['key']]
        for segment in segments[1:]:
            # Merge adjacent segments of same label
            label = mnemonics[segment['key']]
            if label == prev_label and \
                    segment['start'] == prev_segment['end']:
                try:
                    strands_match = (segment['strand'] == prev_segment['strand'])
                except IndexError:
                    strands_match = True

                if strands_match:
                    prev_segment['end'] = segment['end']
                    continue

            print_segment(prev_segment, prev_label)
            prev_segment, prev_label = segment, label
        print_segment(prev_segment, prev_label)

    if out is not sys.stdout:
        out.close()

def parse_args(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTATION MNEMONICFILE"
    parser = OptionParser(usage=usage, version=__version__,
                          description=__doc__.strip())
    add_common_options(parser, ['quiet'])
    parser.add_option("-o", "--outfile", dest="outfile",
                      default=None, metavar="FILE",
                      help="Save relabeled bed file to FILE instead of"
                      " printing to stdout (default)")

    options, args = parser.parse_args(args)

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)

    for arg in args:
        if not os.path.isfile(arg):
            die("Could not find file: %s" % arg)

    kwargs = dict(options.__dict__)
    relabel(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
