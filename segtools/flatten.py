#!/bin/env python

"""
Combine segments from all SEGMENTATION files, labeling with unique labels for
each combination of overlapping that occurs.
Outputs a bed file on
stdout (-o to change), and generates a file in the current directory
that describes the generated labels (use -m to change).
"""

# Author: Orion Buske <stasis@uw.edu>
# Date:   02.02.2010


from __future__ import division, with_statement

import os
import sys

from numpy import concatenate

from . import log, Segmentation, die, add_common_options
from .version import __version__

DEFAULT_HELPFILE = "flat.mnemonics"


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
    """Merge together segments from different segmentations

    segmentations: a dict mapping file names to Segmentation objects.
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

def filter_segments(labels, segments, filter=None):
    # Do nothing if no filtering
    if filter is None or filter == 0:
        return labels, segments

    if (not isinstance(filter, float) or filter < 0 or filter > 1):
        raise ValueError("Invalid value of filter: %s" % str(filter))

    # Find span of labels and whole segmentation
    print >>sys.stderr, "Calculating span of labels and segmentation...",
    n_seg = 0
    n_keys = dict([(key, 0) for key in labels.iterkeys()])
    for segment in segments:
        key = segment[3]
        length = segment[2] - segment[1]
        n_keys[key] += length
        n_seg += length
    print >>sys.stderr, "done"

    # Find segment labels that need to be filtered
    filtered_labels = dict(labels)  # Copy of original labels
    thresh = n_seg * filter
    for key, n_key in n_keys.iteritems():
        if n_key < thresh:
            del filtered_labels[key]

    # Remove segments not in filtered segment labels
    print >>sys.stderr, "Filtering segment labels below threshold...",
    filtered_segments = []
    for segment in segments:
        key = segment[3]
        if key in filtered_labels:
            filtered_segments.append(segment)
    print >>sys.stderr, "done"

    return filtered_labels, filtered_segments

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

def print_readme(labels, filename=DEFAULT_HELPFILE):
    if os.path.isfile(filename):
        log("Warning: overwriting file: %s" % filename)

    with open(filename, "w") as ofp:
        ofp.write("%s\n" % "\t".join(["old", "new", "description"]))
        for key, label in labels.iteritems():
            ofp.write("%s\n" % "\t".join([str(key), str(key), label]))

def flatten(files, outfile=None, helpfile=DEFAULT_HELPFILE,
            filter=None, verbose=True):
    segmentations = {}
    for file in files:
        assert os.path.isfile(file)
        nice_filename = os.path.basename(file)
        segmentations[nice_filename] = Segmentation(file, verbose=verbose)

    labels, segments = merge_segments(segmentations)
    labels, segments = filter_segments(labels, segments, filter=filter)
    print_bed(segments, filename=outfile)
    print_readme(labels, filename=helpfile)

def parse_args(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTATION..."
    parser = OptionParser(usage=usage, version=__version__,
                          description=__doc__.strip())

    add_common_options(parser, ['quiet'])
    parser.add_option("-m", "--mnemonic-file", dest="helpfile",
                      default=DEFAULT_HELPFILE, metavar="FILE",
                      help="Save mapping information to FILE instead of"
                      " %default (default). This file complies with the"
                      " mnemonic file format.")
    parser.add_option("-o", "--outfile", dest="outfile",
                      default=None, metavar="FILE",
                      help="Save flattened bed file to FILE instead of"
                      " printing to stdout (default)")
    parser.add_option("-f", "--filter", dest="filter",
                      default=None, type="float", metavar="F",
                      help="Don't output new segment labels (and corresponding"
                      " segments) that span less than F*N bases, where"
                      " N is the number of bases covered by the new"
                      " segmentation. This can be used to remove"
                      " extremely uncommon labels (e.g. F = 0.01)"
                      " that are the more likely to be spurious."
                      " Filtering is off by default.")

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
              "verbose": options.verbose,
              "filter": options.filter,
              "outfile": options.outfile}
    flatten(args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
