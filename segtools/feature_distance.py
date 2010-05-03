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

from numpy import array, NaN

from . import Segmentation
from .common import FeatureScanner, load_features

EXT_BED = "bed"
EXT_GFF = "gff"
EXT_GTF = "gtf"
EXT_GZ = "gz"

STANDARD_FIELDS = ["chrom", "start", "end", "name"]
DELIM = "\t"

def get_file_ext(filename):
    # Ignore g'zipping
    head, tail = os.path.splitext(os.path.basename(filename))
    if tail == os.extsep + EXT_GZ:
        head, tail = os.path.splitext(head)

    if tail:
        return tail[1:]  # Drop leading "."
    else:
        return None  # Didn't find ext

def load_data(filename, verbose=True):
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
    """Print distance line, allowing individual distances to be lists

    If a single distance is a list, multiple lines are printed, one with
    each distance.

    """
    # If any distances are lists, recurse
    for i, distance in enumerate(distances):
        if isinstance(distance, list):
            for sub_distance in distance:
                new_distances = list(distances)
                # Print line for each single value
                new_distances[i] = sub_distance
                print_line(labels, chrom, segment, new_distances, extra=extra)
            return

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

def indices_to_bools(indices, length):
    indices = set([int(index) for index in indices])
    bools = [False] * length
    for index in indices:
        if index < 0 or index >= length:
            raise ValueError(index)
        else:
            bools[index] = True

    return bools

def feature_distance(segment_file, feature_files, verbose=True,
                     correct_strands=[], print_nearest=[], all_overlapping=[]):
    # Lists are of integer file indices (0-indexed), for which files that
    # property applies to.
    try:
        correct_strand = indices_to_bools(correct_strands, len(feature_files))
    except ValueError, e:
            die("File index for strand correction: %s is invalid" % e)

    try:
        print_nearest = indices_to_bools(print_nearest, len(feature_files))
    except ValueError, e:
        die("File index for print nearest: %s is invalid" % e)

    try:
        all_overlapping = indices_to_bools(all_overlapping, len(feature_files))
    except ValueError, e:
        die("File index for all overlapping: %s is invalid" % e)

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
            for file_i, feature_scanner in enumerate(feature_scanners):
                nearest = feature_scanner.nearest(segment)
                distance = NaN
                if nearest is not None:
                    if isinstance(nearest, list):
                        # If we are only printing one, use the first regardless
                        if not all_overlapping[file_i]:
                            nearest = nearest[0:1]
                    else:
                        # Make sure nearest is a list, even if only one element
                        nearest = [nearest]

                    # Save multiple distances as a list of them
                    distance = []
                    for feature in nearest:
                        one_dist = FeatureScanner.distance(segment, feature)
                        if correct_strand[file_i]:
                            try:
                                strand = feature['strand']
                                if strand == ".":
                                    strand = None
                            except IndexError:
                                strand = None

                            if strand is None:
                                die("Trying to use strand information from:"
                                    " %r, but it was not found" %
                                    os.path.basename(feature_files[file_i]))

                            if (strand == "+" and \
                                    feature['start'] < segment['start']) or \
                                    (strand == "-" and \
                                         feature['end'] > segment['end']):
                                one_dist = -one_dist

                        if print_nearest[file_i]:
                            name = feature_objs[file_i].labels[feature['key']]
                            one_dist = "%s %s" % (one_dist, name)

                        distance.append(one_dist)

                distances.append(distance)

            print_line(labels, chrom, segment, distances)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] SEGMENTFILE FEATUREFILE..."
    version = "%%prog %s" % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=__doc__.strip())

    parser.add_option("-q", "--quiet", dest="verbose",
                      default=True, action="store_false",
                      help="Do not print diagnostic messages.")
    parser.add_option("-s", "--strand-correct", dest="correct_strands",
                      default=[], action="append", metavar="N", type="int",
                      help="Strand correct features from the Nth feature file"
                      " (where N=1 is the first file)."
                      " If the feature list contain strand information,"
                      " the sign of the distance value is used to portray the"
                      " relationship between the segment and the nearest"
                      " stranded feature. A positive distance value indicates"
                      " that the segment is nearest the 5' end of the feature,"
                      " whereas a negative value indicates the segment is"
                      " nearest the 3' end of the feature.")
    parser.add_option("-p", "--print-nearest", dest="print_nearest",
                      default=[], action="append", metavar="N", type="int",
                      help="The name of the nearest feature will be printed"
                      " after each distance (with a space separating the"
                      " two) for features from the Nth feature file (where"
                      " N=1 is the first file). If multiple features"
                      " are equally near (or overlap), it is undefined"
                      " which will be printed")
    parser.add_option("-a", "--all-overlapping", dest="all_overlapping",
                      default=[], action="append", metavar="N", type="int",
                      help="If multiple features in the Nth file"
                      " overlap a segment, a separate line"
                      " is printed for each overlapping segment-feature"
                      " pair (where N=1 is the first file). This is most"
                      " interesting with --print-nearest=N. Otherwise,"
                      " the first overlapping segment will be used for"
                      " printing.")

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)

    segment_file = args[0]
    feature_files = args[1:]

    # Convert file indices to 0-indexed
    correct_strands = [int(val) - 1 for val in options.correct_strands]
    print_nearest = [int(val) - 1 for val in options.print_nearest]
    all_overlapping = [int(val) - 1 for val in options.all_overlapping]

    kwargs = {"verbose": options.verbose,
              "correct_strands": correct_strands,
              "print_nearest": print_nearest,
              "all_overlapping": all_overlapping}
    feature_distance(segment_file, feature_files, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
