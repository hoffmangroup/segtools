#!/bin/env python

"""
Combine segments from all segment files, labeling with unique labels for
each combination of overlapping that occurs. Outputs a bed file on
stdout, and generates a README file in the current directory (unless
otherwise specified) that describes the generated labels.
"""

# Author: Orion Buske <stasis@uw.edu>
# Date:   02.02.2010


from __future__ import division, with_statement

__version__ = "$Revision: "

import os
import sys

from numpy import concatenate

from . import Segmentation
from .common import die


DEFAULT_HELPFILE = "README"


class IntervalSegment(object):
    def __init__(self, segment=None, start=None, end=None, keys=[]):
        if segment is None:
            self.start = start
            self.end = end
            self.keys = keys
        else:
            self.start = segment['start']
            self.end = segment['end']
            self.keys = [segment['key']]

    def signature(self):
        return tuple(sorted(set(self.keys)))

    def __str__(self):
        return "[%r, %r, %s]" % (self.start, self.end, self.keys)

    def __repr__(self):
        return "<IntervalSegment %s>" % str(self)

class Interval(object):
    def __init__(self):
        self.start = None
        self.end = None
        self._segments = []  # Not kept sorted

    def add(self, segment):
        """Add a segment to the interval"""
        if not self._segments:
            self.start = segment['start']
            self.end = segment['end']
            self._segments.append(IntervalSegment(segment=segment))
        else:
            segment_start = segment['start']
            segment_end = segment['end']
            segment_key = segment['key']
            assert segment_start < self.end

            # Potentially partition existing segments along boundaries
            segments_to_add = []
            for i, i_segment in enumerate(self._segments):
                i_start = i_segment.start
                i_end = i_segment.end
                i_keys = i_segment.keys
                if i_start < segment_end and \
                        i_end > segment_start:

                    # Maybe add new segment before overlap
                    if i_start < segment_start:
                        new_segment = IntervalSegment(start=i_start,
                                                      end=segment_start,
                                                      keys=i_keys)
                        segments_to_add.append(new_segment)

                    # Change segment to overlap boundary
                    i_segment.start = max(segment_start, i_start)
                    i_segment.end = min(segment_end, i_end)
                    i_segment.keys = i_keys + [segment_key]

                    # Maybe add new segment after overlap
                    if i_end > segment_end:
                        new_segment = IntervalSegment(start=segment_end,
                                                      end=i_end,
                                                      keys=i_keys)
                        segments_to_add.append(new_segment)

            # Maybe add one new segment at the end
            if segment_end > self.end:
                new_segment = IntervalSegment(start=self.end,
                                              end=segment_end,
                                              keys=[segment_key])
                segments_to_add.append(new_segment)
                self.end = segment_end  # Extend interval end

            # Add all new segments to the end of the interval
            self._segments.extend(segments_to_add)


    def flush_to(self, pos):
        """Return segments from the start of the interval to the position"""
        if not self._segments or pos <= self.start:
            return []
        else:
            segments = []
            i = 0
            while i < len(self._segments):
                segment = self._segments[i]
                if segment.end <= pos:
                    # Segment to be completely flushed
                    segments.append(segment)
                    del self._segments[i]
                else:
                    i += 1

            self.start = pos
            return segments

    def flush(self):
        """Flush all remaining segments"""
        self.start = self.end
        segments = self._segments
        self._segments = []
        return segments

    def __repr__(self):
        return "<Interval: %r, %r, %r>" % \
            (self.start, self.end, self._segments)

def join_chrom_segments(chrom, files, segmentations, label_offsets):
    # Combine segments from all segmentations
    segments = None
    for file in files:
        segmentation = segmentations[file]
        label_offset = label_offsets[file]
        try:
            cur_segments = segmentation.chromosomes[chrom]

            # Modify label keys to keep different file segmentations apart
            cur_segments['key'] += label_offset

            # Add current segment array to segments
            if segments is None:
                segments = cur_segments
            else:
                segments = concatenate((segments, cur_segments), axis=0)
        except KeyError:
            pass

    segments.sort()
    return segments

def flatten_segments(segments):
    """Flatten overlapping segments into IntervalSegments"""

    new_segments = []
    if len(segments) == 0:
        return new_segments

    interval = Interval()
    for segment in segments:
        # Process segments up to segment_start
        flushed_segments = interval.flush_to(segment['start'])
        new_segments.extend(flushed_segments)

        # Add new segment to interval
        interval.add(segment)

    # Flush remaining segments from interval
    new_segments.extend(interval.flush())

    return new_segments

def merge_segments(segmentations):
    """Merges together segments from different segmentations

    segmentations: a dict mapping bed file names to Segmentation objects.
    """

    files = sorted(segmentations.keys())  # Constant file order

    chroms = set()
    label_offsets = {}  # file -> label_offset
    shifted_labels = {}  # shifted_key -> label_string
    for file in files:
        segmentation = segmentations[file]
        chroms.update(segmentation.chromosomes.keys())
        labels = segmentation.labels

        # Shift labels of file to be unique
        label_offset = len(shifted_labels)
        label_offsets[file] = label_offset
        for key, label in labels.iteritems():
            shifted_labels[key + label_offset] = "%s:%s" % (file, label)

    signature_labels = {}  # (shifted_key, ...) -> flat_key
    flat_labels = {}  # flat_key -> label_string
    new_segments = []
    for chrom in chroms:
        segments = join_chrom_segments(chrom, files, segmentations,
                                       label_offsets)
        flat_segments = flatten_segments(segments)

        for flat_segment in flat_segments:
            signature = flat_segment.signature()
            try:
                key = signature_labels[signature]
            except KeyError:
                key = len(signature_labels)
                key_str = "[%s]" % ", ".join([shifted_labels[shifted_key]
                                              for shifted_key in signature])
                signature_labels[signature] = key
                flat_labels[key] = key_str

            new_segment = (chrom, flat_segment.start, flat_segment.end, key)
            new_segments.append(new_segment)

    return flat_labels, new_segments

def print_bed(segments, filename=None):
    """Print segments in BED format to file (or stdout if None)"""
    if filename is None:
        outfile = sys.stdout
    else:
        if os.path.isfile(filename):
            print >>sys.stderr, "Warning: overwriting file: %s" % filename

        outfile = open(filename, "w")

    for segment in segments:
        print >>outfile, "\t".join([str(val) for val in segment])

    if outfile is not sys.stdout:
        outfile.close()

def print_readme(labels, filename=DEFAULT_HELPFILE):
    if os.path.isfile(filename):
        print >>sys.stderr, "Warning: overwriting file: %s" % filename

    with open(filename, "w") as ofp:
        print >>ofp, "# new_label\toverlap(file:old_label, ...)"
        for key, label in labels.iteritems():
            print >>ofp, "\t".join([str(key), label])

def flatten_bed(bedfiles, outfile=None, helpfile=DEFAULT_HELPFILE):
    segmentations = {}
    for bedfile in bedfiles:
        assert os.path.isfile(bedfile)
        nice_filename = os.path.basename(bedfile)
        segmentations[nice_filename] = Segmentation(bedfile)

    labels, segments = merge_segments(segmentations)
    print_bed(segments, filename=outfile)
    print_readme(labels, filename=helpfile)

def parse_args(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTS..."
    parser = OptionParser(usage=usage, version=__version__,
                          description=__doc__.strip())
    parser.add_option("-H", "--helpfile", dest="helpfile",
                      default=DEFAULT_HELPFILE, metavar="FILE",
                      help="Save mapping information to FILE instead of"
                      " %default (default).")
    parser.add_option("-o", "--outfile", dest="outfile",
                      default=None, metavar="FILE",
                      help="Save flattened bed file to FILE instead of"
                      " printing to stdout (default)")
    options, args = parser.parse_args(args)

    if len(args) == 0:
        parser.error("Inappropriate number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)

    for arg in args:
        if not os.path.isfile(arg):
            die("Could not find file: %s" % arg)

    kwargs = {"helpfile": options.helpfile,
              "outfile": options.outfile}
    flatten_bed(args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
