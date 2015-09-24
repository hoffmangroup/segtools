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

import os
import sys

from . import log, Segmentation, die, add_common_options
from .common import map_mnemonics, get_ordered_labels
from .version import __version__

class StrandUnknown(str):
    """sentinel object that converts to '.'"""

    def __str__(self):
        return "."

SCORE_DEFAULT = 0
STRAND_UNKNOWN = StrandUnknown()

def get_strand(segment):
    try:
        return segment['strand']
    # IndexError, not KeyError, because segment is a NumPy struct
    except IndexError:
        return STRAND_UNKNOWN

def relabel(segfilename, mnemonicfilename, outfile=None, verbose=True):
    assert os.path.isfile(segfilename)
    segmentation = Segmentation(segfilename, verbose=verbose)

    labels = segmentation.labels

    mnemonics = map_mnemonics(labels, mnemonicfilename)
    ordered_labels, mnemonics = get_ordered_labels(labels, mnemonics)

    try:
        colors = map_mnemonics(labels, mnemonicfilename, "color")
        ordered_colors, colors = get_ordered_labels(labels, colors)
    except KeyError:
        colors = None

    log("Found labels:\n%s" % labels)
    log("Found mnemonics:\n%s" % mnemonics)
    log("Found colors:\n%s" % colors)

    if outfile is None:
        out = sys.stdout
    else:
        out = open(outfile, "w")

    try:
        for chrom, segments in segmentation.chromosomes.iteritems():
            def print_segment(segment, label):
                start = segment['start']
                end = segment['end']
                old = segment['key']
                tokens = [chrom, start, end, mnemonics[old]]

                strand = get_strand(segment)
                if strand != STRAND_UNKNOWN or colors:
                    tokens.extend([SCORE_DEFAULT, str(strand)])

                if colors:
                    tokens.extend([start, end, colors[old]])
                print >>out, "\t".join(map(str, tokens))

            prev_segment = segments[0]
            prev_label = mnemonics[prev_segment['key']]
            for segment in segments[1:]:
                # Merge adjacent segments of same label
                label = mnemonics[segment['key']]
                if (label == prev_label and
                    segment['start'] == prev_segment['end']):
                    try:
                        strands_match = \
                            (segment['strand'] == prev_segment['strand'])
                    except IndexError:
                        strands_match = True

                    if strands_match:
                        prev_segment['end'] = segment['end']
                        continue

                print_segment(prev_segment, prev_label)
                prev_segment, prev_label = segment, label
            print_segment(prev_segment, prev_label)
    # XXX: this would be cleaner as a context manager
    finally:
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
