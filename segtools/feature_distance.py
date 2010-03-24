#!/usr/bin/env python

"""
Given a list of segments and lists of features, prints the distance to
the nearest feature in each list (zero if overlap), for each segment.
Distance is the difference between the nearest bases of the segment and the,
features so features with one base pair between them are a distance of two
apart.

Segments and features should be specified in [gzip'd] files in either
BED3+ or GFF format.
"""

from __future__ import division, with_statement

__version__ = "$Revision$"


import os
import sys

from numpy import array, NaN, ndarray
from warnings import warn

from . import Segmentation
from .common import load_features

EXT_BED = "bed"
EXT_GFF = "gff"
EXT_GTF = "gtf"
EXT_GZ = "gz"

STANDARD_FIELDS = ["chrom", "start", "end", "name"]
DELIM = "\t"

class FeatureScanner(object):
    def __init__(self, features):
        if isinstance(features, ndarray):
            self._iter = iter(features)
            self._has_next = features.shape[0] > 0
        else:
            warn("Feature list not of understood type. Treated as empty")
            self._has_next = False

        self._prev_best = None
        self._next_best = None
        self._prev_start = None

    @staticmethod
    def distance(feature1, feature2):
        """Return the distance between two features (zero if overlap)"""
        if feature1['end'] <= feature2['start']:
            return feature2['start'] - feature1['end'] + 1
        elif feature2['end'] <= feature1['start']:
            return feature1['start'] - feature2['end'] + 1
        else:  # Features overlap
            return 0

    def should_iter(self, segment):
        if not self._has_next:
            return False  # No more features to read
        elif self._prev_best is None:
            return True  # Haven't read two features yet
        else:
            return self._next_best['start'] < segment['end']

    def next_segment_nearest(self, segment):
        """Return the nearest feature for the next segment

        Should be called with in-order segments

        """
        while self.should_iter(segment):
            # Update previous nearest
            if self._prev_best is None or self._next_best is None or \
                self._next_best['end'] >= self._prev_best['end']:
                self._prev_best = self._next_best

            # Advance next nearest
            try:
                self._next_best = self._iter.next()

                # Ensure in-order segments
                next_start = self._next_best['start']
                if self._prev_start is not None and \
                        self._prev_start > next_start:
                    raise ValueError("Segments provided out of order."
                                     " Last segment ended at %r, but current"
                                     " one ends at %r." % \
                                         (self._prev_start, next_start))
                else:
                    self._prev_start = next_start
            except StopIteration:
                self._has_next = False

        if self._prev_best is None or self._next_best is None:
            if self._prev_best is None and self._next_best is None:
                # There are no features to compare against
                return None
            elif self._prev_best is None:
                return self._next_best
            else:  #self._next_best is None
                return self._prev_best
        else:
            # Return best of both pointers
            if self.distance(self._prev_best, segment) <= \
                    self.distance(self._next_best, segment):
                return self._prev_best
            else:
                return self._next_best

def get_file_ext(filename):
    # Ignore g'zipping
    head, tail = os.path.splitext(os.path.basename(filename))
    if tail == os.extsep + EXT_GZ:
        head, tail = os.path.splitext(head)

    if tail:
        return tail[1:]  # Drop leading "."
    else:
        return None  # Didn't find ext

def load_data(filename, verbose=False):
    ext = get_file_ext(filename)
    if ext == EXT_BED:
        return Segmentation(filename, verbose=verbose)
    elif ext == EXT_GFF or ext == EXT_GTF:
        return load_features(filename, verbose=verbose)
    else:
        raise NotImplementedError("Unexpected file type: %s" % ext)

def print_header_line(feature_files):
    file_fields = [os.path.basename(filename) for filename in feature_files]
    fields = STANDARD_FIELDS + file_fields
    print DELIM.join(fields)

def print_line(labels, chrom, segment, distances, extra=[]):
    distance_strs = [str(distance) for distance in distances]
    segment_strs = [chrom,
                    str(segment['start']),
                    str(segment['end']),
                    labels[segment['key']]]
    extra_strs = [str(extra_field) for extra_field in extra]
    fields = segment_strs + distance_strs + extra_strs
    print DELIM.join(fields)

def die(msg):
    raise Exception("ERROR: %s" % msg)

def feature_distance(segment_file, feature_files, verbose=False,
                     correct_strands=[], print_nearest=False):
    correct_strands = set([int(col) for col in correct_strands])
    correct_strand = [False] * len(feature_files)
    for col in correct_strands:
        if col < 0 or col >= len(feature_files):
            die("Strand correction column: %d is invalid" % col)
        else:
            correct_strand[col] = True

    # Load segment file
    segment_obj = load_data(segment_file, verbose=verbose)
    labels = segment_obj.labels
    segment_data = segment_obj.chromosomes

    # Load feature files
    feature_objs = [load_data(feature_file, verbose=verbose)
                    for feature_file in feature_files]
    feature_datas = [feature_obj.chromosomes for feature_obj in feature_objs]

    print_header_line(feature_files)
    for chrom in segment_data:
        if verbose:
            print >>sys.stderr, "%s" % chrom

        segments = segment_data[chrom]
        features_list = []
        for feature_data in feature_datas:
            try:
                features = feature_data[chrom]
            except KeyError:
                features = array(())
            features_list.append(features)

        feature_scanners = [FeatureScanner(features)
                            for features in features_list]
        for segment in segments:
            distances = []
            for features_i, feature_scanner in enumerate(feature_scanners):
                feature = feature_scanner.next_segment_nearest(segment)
                distance = NaN
                if feature is not None:
                    distance = FeatureScanner.distance(segment, feature)
                    col_correct_strand = correct_strand[features_i]
                    if col_correct_strand:
                        try:
                            strand = feature['strand']
                            if strand == ".":
                                strand = None
                        except IndexError:
                            strand = None

                        if strand is None:
                            die("Trying to use strand information from: %r,"
                                " but it was not found" % os.path.basename(
                                    feature_files[features_i]))

                        if (strand == "+" and \
                                feature['start'] < segment['start']) or \
                                (strand == "-" and \
                                     feature['end'] > segment['end']):
                            distance = -distance

                    if print_nearest:
                        name = feature_objs[features_i].labels[feature['key']]
                        distance = "%s %s" % (distance, name)

                distances.append(distance)

            print_line(labels, chrom, segment, distances)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTFILE FEATUREFILE..."
    version = "%%prog %s" % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=__doc__.strip())

    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print diagnostic messages")
    parser.add_option("-s", "--strand-correct", dest="correct_strands",
                      default=[], action="append", metavar="N",
                      help="Strand correct features from the Nth feature file"
                      " (where N=0 is the first file)."
                      " If the feature list contain strand information,"
                      " the sign of the distance value is used to portray the"
                      " relationship between the segment and the nearest"
                      " stranded feature. A positive distance value indicates"
                      " that the segment is nearest the 5' end of the feature,"
                      " whereas a negative value indicates the segment is"
                      " nearest the 3' end of the feature.")
    parser.add_option("-p", "--print-nearest", dest="print_nearest",
                      action="store_true", default=False,
                      help="The name of the nearest feature will be printed"
                      " after each distance (with a space separating the"
                      " two).")

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)

    segment_file = args[0]
    feature_files = args[1:]
    kwargs = {"verbose": options.verbose,
              "correct_strands": options.correct_strands,
              "print_nearest": options.print_nearest}
    feature_distance(segment_file, feature_files, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
