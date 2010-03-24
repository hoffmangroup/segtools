#!/usr/bin/env python
from __future__ import division, with_statement

"""
validate.py:

Provides command-line and package entry points for analyzing the signal
distribution over tracks and labels.

"""

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE="signal_distribution"

__version__ = "$Revision$"


import os
import sys

from collections import defaultdict
from genomedata import Genome
from itertools import repeat
from functools import partial
from numpy import array, ceil, floor, fromiter, histogram, isfinite, NAN, nanmax, nanmin, nansum, NINF, PINF
from rpy2.robjects import r, numpy2ri

from . import Segmentation
from .common import die, image_saver, iter_segments_continuous, make_tabfilename, map_mnemonics, r_source, setup_directory, tab_reader, tab_saver
from .html import save_html_div
from .mnemonics import create_mnemonic_file

FIELDNAMES = ["label", "trackname", "lower_edge", "count"]
NAMEBASE = "%s" % MODULE

FIELDNAMES_STATS = ["label", "trackname", "mean", "sd"]
NAMEBASE_STATS = os.extsep.join([NAMEBASE, "stats"])

HTML_TITLE = "Signal value distribution"
HTML_TEMPLATE_FILENAME = "signal_div.tmpl"

PNG_WIDTH_PER_SUBPLOT = 200
PNG_HEIGHT_PER_SUBPLOT = 150
PNG_BASE_HEIGHT = 200
PNG_BASE_WIDTH = 100

STATS_PNG_SIZE = 1000

NBINS = 100

class SignalHistogram(object):
    def __init__(self, histogram=None, **fields):
        self._data = histogram
        self._metadata = fields

    @property
    def data(self):
        return self._data

    @property
    def metadata(self):
        return self._metadata

    @staticmethod
    def calculate(genome, segmentation, nbins=None, chroms=None,
                  group_labels=False, calc_ranges=False,
                  value_range=(None, None), quick=False, verbose=True):
        """Computes a set of histograms, one for each track-label pair.

        if nbins is:
        - None: There will be a bin for every possible value
        - > 0: There will be this many bins

        if chroms is:
        - a sequence of chromosome names: only those chromosomes are processed
        - False, [], None: all chromosomes are processed

        if group_labels is:
        - True: Group together all labels and only generate one histogram
          for the track over all labels
        - False: Generate a histogram for every track-label pair

        if calc_ranges is:
        - False: use the precomputed ranges stored for the entire genome file
        - True: calculate the segmentation ranges across the entire input

        value_range is a (min, max) tuple:
          if either value is None, it is ignored
          if min/max is not None, it is used as a limit for all binning
          if both are specified, calc_ranges is ignored

        if quick is:
        - True: the histogram is calculated for only the first chromosome

        returns a dict

        if group_labels:
          key: trackname
          val: [hist, edges]
        else:
          key: label (not label_key)
          val: dict
               key: trackname
               val: [hist, edges]
        """
        assert genome is not None
        assert segmentation is not None

        labels = segmentation.labels

        label_bins = xrange(0, len(labels) + 1)  # [0, 1, ..., n, n + 1]

        tracks = genome.tracknames_continuous
        max_bins = ceil(genome.maxs).astype(int)
        min_bins = floor(genome.mins).astype(int)
        # A dict from tracks to a range tuple
        track_ranges = dict(zip(tracks, zip(min_bins, max_bins)))

        if len(tracks) == 0:
            die("Trying to calculate histogram for no tracks")

        # key: trackname
        # val: track_hist_func: a function that generates a histogram from data
        #                       with uniform bins within a track
        # For every trackname, there's a function that generates histograms
        track_hist_funcs = dict((trackname,
                                 partial(histogram,
                                         bins=xrange(min_bin, max_bin + 1),
                                         new=True))
                                for trackname, (min_bin, max_bin) in
                                track_ranges.iteritems())
        # Trackname is one the track names, track_range is (min,max)
        # functools.partial(func[, *args][, **keywords])
        #   Returns a new partial object which when called will behave like
        #   func called with args and keywords.

        #if group_labels:  # Pretend there is just one "all" label
        #    labels = dict([(-1, "all")])

        # Result dictionary: see function docstring for specification
        # Now for every trackname and for every label
        #   there's a histogram array of the right boundary and frequency
        res = dict((trackname, dict((label,
                                     list(track_hist_func(array([]))))
                                    for label in labels.itervalues()))
                   for trackname, track_hist_func in \
                       track_hist_funcs.iteritems())

        nseg_dps = 0  # Number of non-NaN data points in segmentation tracks
        if verbose:
            print >>sys.stderr, "Generating signal distribution histograms"

        with genome:
            if chroms:
                chromosomes = [genome[chrom] for chrom in chroms]
            else:
                chromosomes = genome

            for chromosome in chromosomes:
                # Iterate through supercontigs and segments together
                for segment, continuous_seg in \
                        iter_segments_continuous(chromosome, segmentation,
                                                 verbose=verbose):
                    seg_label = labels[segment['key']]
                    # Iterate through each track
                    for trackname, track_hist_func in \
                            track_hist_funcs.iteritems():
                        # col_index is the index of the trackname
                        #   (if 5 tracks, it can be 0, 1, 2, 3 or 4)
                        col_index = chromosome.index_continuous(trackname)

                        # continuous_col is the "intensity"
                        #   (continuous number) from the data track
                        # len(supercontig_map) = len(continuous_col)
                        #   because they are for the same segment
                        continuous_col = continuous_seg[:, col_index]

                        # Remove the NaN's, otherwise numpy.histogram
                        #   gets confused and misses some of the data
                        cur_col_nonan = remove_nans(continuous_col)

                        # Keep track of number of real data values processed
                        if trackname in segmentation.tracks:
                            nseg_dps += cur_col_nonan.shape[0]

                        # Edges is a list with the left edge of every bin,
                        #   hist is the frequency for every bin
                        # This calls numpy.histogram
                        hist, edges = track_hist_func(cur_col_nonan)

                        cur_res = res[trackname][seg_label]
                        cur_res[0] += hist

                if quick: break  # 1 chromosome

        if verbose:
            print >>sys.stderr, "Segmentation tracks: %s\n" % \
                str(segmentation.tracks)
            print >>sys.stderr, ("Read %s non-NaN values from"
                                 " segmentation tracks" % nseg_dps)

        return SignalHistogram(res, nseg_dps=nseg_dps)

    @staticmethod
    def read(dirpath, basename=NAMEBASE, **kwargs):
        """Read a histogram from an output directory"""
        with tab_reader(dirpath, basename, **kwargs) as (reader, metadata):
            histogram = {}
            try:
                nseg_dps = int(metadata["nseg_dps"])
            except KeyError:
                nseg_dps = 0

            for row in reader:
                type = row[0]
                if type == "track":
                    trackname = row[1]
                    edges = fromiter(row[2:], float)
                    histogram[trackname] = {}
                elif type == "label":
                    label = row[1]
                    hist = fromiter(row[2:], int)
                    histogram[trackname][label] = (hist, edges)
                else:
                    raise IOError("Unexpected row label: %s in input"
                                  " directory: %s" % (type, dirpath))

        return SignalHistogram(histogram, nseg_dps=nseg_dps)

    def save(self, dirpath, basename=NAMEBASE, **kwargs):
        """Saves the histogram data to a tab file"""
        with tab_saver(dirpath, basename, metadata=self.metadata,
                       **kwargs) as saver:
            for trackname in sorted(self.data.keys()):
                label_hists = self.data[trackname]
                first = True
                for label in sorted(label_hists.keys()):
                    (hist, edges) = label_hists[label]
                    if first:
                        nice_edges = edges.tolist()
                        saver.writerow(["track", trackname] + nice_edges)
                        first = False

                    nice_hist = hist.astype(int).tolist()
                    saver.writerow(["label", label] + nice_hist)

    def add(self, o_histogram):
        """Add data from a second histogram"""
        assert isinstance(o_histogram, SignalHistogram)
        if self.data is None:
            self._data = o_histogram.data
        else:
            for trackname, label_hists in self.data.iteritems():
                o_label_hists = o_histogram.data[trackname]
                for label, (hist, edges) in label_hists.iteritems():
                    o_hist, o_edges = o_label_hists[label]
                    assert (edges == o_edges).all()
                    hist += o_hist

def calc_stats(histogram):
    """Calculate track statistics (mean, sd) for each label

    histogram: a SignalHistogram

    Values are approximated from binned distributions
    - mean is a lower-bound on the actual mean

    Returns a dict: label_key -> dict( trackname -> {"mean", "sd", ...} )

    """
    stats = defaultdict(partial(defaultdict, dict))
    for trackname, label_hists in histogram.data.iteritems():
        for label, (hist, edges) in label_hists.iteritems():
            cur_stat = stats[label][trackname]
            n = nansum(hist)
            if n == 0:
                mean = NAN
                sd = NAN
            else:
                mean = (hist * edges[:-1]).sum() / n
                sd = (hist * (edges[:1] - mean)**2).sum() / (n - 1)

            cur_stat["sd"] = sd
            cur_stat["mean"] = mean

    return stats


## Saves the track stats to a tab file
def save_stats_tab(stats, dirpath, clobber=False,
                   basename=NAMEBASE_STATS, fieldnames=FIELDNAMES_STATS):
    with tab_saver(dirpath, basename, fieldnames, clobber=clobber) as saver:
        for label, label_stats in stats.iteritems():
            for trackname, track_stats in label_stats.iteritems():
                mean = track_stats["mean"]
                sd = track_stats["sd"]
                saver.writerow(locals())

def start_R():
    r_source("common.R")
    r_source("signal.R")

def constant_factory(val):
    return repeat(val).next

def remove_nans(numbers):
    return numbers[isfinite(numbers)]

## Returns a dict of trackname -> tuples: (min, max), one for each trackname
def seg_min_max(chromosome, segmentation, verbose=False):
    #print >>sys.stderr, "Looking for min & max in chromosome ", \
    #     chromosome.name, " for track ", track_index
    #print >>sys.stderr, "segmentation: ", segmentation.chromosomes

    tracknames = chromosome.tracknames_continuous
    limits = dict([(trackname, (PINF, NINF)) for trackname in tracknames])

    if chromosome.name not in segmentation.chromosomes.keys():
        return limits

    for segment, continuous in \
            iter_segments_continuous(chromosome, segmentation,
                                     verbose=verbose):
        if continuous.shape[0] > 0:
            for trackname in tracknames:
                col_index = chromosome.index_continuous(trackname)
                continuous_col = continuous[:, col_index]

                if len(continuous_col) != 0:
                    old_min, old_max = limits[trackname]
                    cur_min = nanmin(continuous_col, axis=0)
                    cur_max = nanmax(continuous_col, axis=0)
                    limits[trackname] = (min(old_min, cur_min),
                                         max(old_max, cur_max))

    return limits

## Loads the ranges of each track from the genomedata object
##   unless a segmentation is specified, in which case the ranges
##   are calculated for that segmentation
## Thus, defaults to loading ranges for all tracks from genomedata
def load_track_ranges(genome, segmentation=None):
    """
    returns a dict
    key: trackname
    val: (min, max) for that trackname
    """
    # start with the most extreme possible values
    res = defaultdict(constant_factory((PINF, NINF)))
    for chromosome in genome:
        if segmentation is not None:
            print >>sys.stderr, "\t%s" % chromosome.name
            res_chrom = seg_min_max(chromosome, segmentation)
            for trackname, limits in res_chrom.iteritems():
                old_min, old_max = res[trackname]
                cur_min, cur_max = limits
                res[trackname] = (min(old_min, cur_min), max(old_max, cur_max))
                #print >>sys.stderr, "Limits for %s: %s" % \
                #     (trackname, res[trackname])
        else:
            tracknames = chromosome.tracknames_continuous
            for trackname in tracknames:
                col_index = chromosome.index_continuous(trackname)
                cur_min = chromosome.mins[col_index]
                cur_max = chromosome.maxs[col_index]

                old_min, old_max = res[trackname]
                res[trackname] = (min(old_min, cur_min), max(old_max, cur_max))

    # Cast out of defaultdict
    return dict(res)

# ## Computes a set of histograms, one for each track-label pair.
# def calc_histogram(genome, segmentation, nbins=None,
#                    group_labels=False, calc_ranges=False,
#                    value_range=(None, None), quick=False):
#     """
#     if nbins is:
#     - None: There will be a bin for every possible value
#     - > 0: There will be this many bins

#     if group_labels is:
#     - True: Group together all labels and only generate one histogram
#       for the track over all labels
#     - False: Generate a histogram for every track-label pair

#     if calc_ranges is:
#     - False: use the precomputed ranges stored for the entire genome file
#     - True: calculate the segmentation ranges across the entire input

#     value_range is a (min, max) tuple:
#       if either value is None, it is ignored
#       if min/max is not None, it is used as a limit for all binning
#       if both are specified, calc_ranges is ignored

#     if quick is:
#     - True: the histogram is calculated for only the first chromosome

#     returns a dict

#     if group_labels:
#       key: trackname
#       val: [hist, edges]
#     else:
#       key: label_key
#       val: dict
#            key: trackname
#            val: [hist, edges]
#     """
#     assert genome is not None
#     assert segmentation is not None

#     if calc_ranges and any(lim is None for lim in value_range):
#         # Calculate ranges over segmentation only
#         print >>sys.stderr, "Calculating track ranges...",
#         track_ranges = load_track_ranges(genome, segmentation)
#         print >>sys.stderr, "done."
#     else:
#         # Lookup ranges in genomedata
#         track_ranges = load_track_ranges(genome)

#     # Override ranges as necessary
#     value_min, value_max = value_range
#     for trackname, track_range in track_ranges.iteritems():
#         track_min, track_max = track_range
#         if value_min is not None:
#             track_min = value_min
#         if value_max is not None:
#             track_max = value_max
#         track_ranges[trackname] = (track_min, track_max)

#     tracks = track_ranges.keys()  # Re-select just the tracks used

#     if len(tracks) == 0:
#         die("Trying to calculate histogram for no tracks")

#     # Calculate bins for each track
#     if nbins is None:  # Have a bin for every possible integer value
#         track_bins = dict([(trackname,
#                             math.ceil(track_range[1]) -
#                             math.floor(track_range[0]))
#                            for trackname, track_range in track_ranges])
#     else:  # Use the given number of bins
#         assert nbins > 0
#         track_bins = dict([(trackname, nbins) for trackname in tracks])

#     # key: trackname
#     # val: histogram_func: a function that generates a histogram from data
#     #                      with uniform bins within a track
#     # For every trackname, there's a function that generates histograms
#     histogram_funcs = dict((trackname, partial(histogram,
#                                                bins=track_bins[trackname],
#                                                range=track_ranges[trackname],
#                                                new=True))
#                            for trackname in tracks)
#     # Trackname is one the track names, track_range is (min,max)
#     # functools.partial(func[, *args][, **keywords])
#     #   Returns a new partial object which when called will behave like func
#     #   called with args and keywords.

#     labels = segmentation.labels
#     if group_labels:  # Pretend there is just one "all" label
#         labels = dict([(-1, "all")])

#     # Return dictionary: see function docstring for specification
#     # Now for every label and for every trackname
#     #   there's a histogram array of the right boundary and frequency
#     res = dict((label_key, dict((trackname,
#                                  list(histogram_func(array([]))))
#                                 for trackname, histogram_func
#                                 in histogram_funcs.iteritems()))
#                for label_key in labels)

#     print >>sys.stderr, "Generating signal distribution histograms"
#     ntracks = len(tracks)
#     nlabels = len(labels)
#     nseg_dps = 0  # Number of non-NaN data points in segmentation tracks

#     with genome:
#         for chromosome in genome:
#             chrom = chromosome.name

#             # Iterate through supercontigs and segments together
#             for segment, continuous_seg in \
#                     iter_segments_continuous(chromosome, segmentation):

#                 if group_labels:
#                     cur_res = res  # All labels are used
#                 else:
#                     seg_label_key = segment['key']
#                     cur_res = res[seg_label_key]

#                 # Iterate through each track
#                 for trackname, histogram_func in histogram_funcs.iteritems():
#                     # col_index is the index of the trackname
#                     #   (if 5 tracks, it can be 0, 1, 2, 3 or 4)
#                     col_index = chromosome.index_continuous(trackname)

#                     # continuous_col is the "intensity"
#                     #   (continuous number) from the data track
#                     # len(supercontig_map) = len(continuous_col)
#                     #   because they are for the same segment
#                     continuous_col = continuous_seg[:, col_index]

#                     # Remove the NaN's, otherwise the numpy.histogram function
#                     #   gets confused and misses some of the data
#                     cur_col_nonan = remove_nans(continuous_col)

#                     # Keep track of number of real data values processed
#                     if trackname in segmentation.tracks:
#                         nseg_dps += cur_col_nonan.shape[0]

#                     # Edges is a list containing the left edge of every bin,
#                     #   hist is the frequency for every bin
#                     # This calls numpy.histogram
#                     hist, edges = histogram_func(cur_col_nonan)

#                     res_trackname = cur_res[trackname]

#                     assert (res_trackname[1] == edges).all()
#                     res_trackname[0] += hist
#             if quick: break  # 1 chromosome

#     print ""
#     print "Segmentation tracks: %s" % str(segmentation.tracks)
#     print "Read %s non-NaN values from segmentation tracks\n" % nseg_dps
#     return res, nseg_dps







# ## Computes a set of histograms, one for each track-label pair.
# def calc_histogram2(genome, segmentation, quick=False, verbose=True,
#                    group_labels=False, calc_ranges=False,
#                    value_range=(None, None), nbins=None,):
#     """
#     if quick is:
#     - True: the histogram is calculated for only the first chromosome

#     returns a dict:
#       key: trackname
#       val: list of [hist, label_edges, value_edges]:
#         hist: numpy.ndarray:
#           rows: labels [0, n]
#           cols: value_bins
#         label_edges: numpy.ndarray of label bins
#         value_edges: numpy.ndarray of value bins
#     """
#     assert genome is not None
#     assert segmentation is not None

#     labels = segmentation.labels

#     label_bins = xrange(0, len(labels) + 1)  # [0, 1, ..., n, n + 1]

#     tracks = genome.tracknames_continuous
#     max_bins = ceil(genome.maxs).astype(int)
#     min_bins = floor(genome.mins).astype(int)
#     # A dict from tracks to a range tuple
#     track_ranges = dict(zip(tracks, zip(min_bins, max_bins)))

#     # A dict of functions, one for each track that,
#     #   given (label_map, track_continuous) generates a
#     #   2d-histogram across labels and continuous values
#     track_hist_funcs = dict([(trackname,
#                               partial(histogram2d,
#                                       bins=[label_bins,
#                                             xrange(min_bin, max_bin + 1)]))
#                              for trackname, (min_bin, max_bin) in
#                              track_ranges.iteritems()])

#     # Return data: see function docstring for specification
#     # Initialize to appropriate size and type by calling functions with
#     #   values that won't be counted
#     res = {}
#     for trackname in tracks:
#         hist, xedges, yedges = track_hist_funcs[trackname]([PINF], [PINF])
#         res[trackname] = [hist.astype(int), xedges.astype(int),
#                           yedges.astype(int)]

#     print >>sys.stderr, "Generating signal distribution histograms"
#     ntracks = len(tracks)
#     nlabels = len(labels)
#     # Number of defined (non-NaN) data points found in segmentation tracks
#     nseg_dps = 0


#     nsupercontigs = 0
#     with genome:
#         for chromosome in genome:
#             for supercontig, segments in \
#                     iter_supercontig_segments(chromosome, segmentation,
#                                               verbose=False):
#                 nsupercontigs += 1

#     print "Created progress bar from %d items" % (nsupercontigs * ntracks)
#     progress = ProgressBar(nsupercontigs * ntracks)

#     with genome:
#         for chromosome in genome:
#             for supercontig, segments in \
#                     iter_supercontig_segments(chromosome, segmentation,
#                                               verbose=False):
#                 continuous = supercontig.continuous
#                 start_i = max(segments['start'][0], supercontig.start)
#                 end_i = min(segments['end'][-1], supercontig.end)
#                 label_map = map_segment_label(segments, (start_i, end_i))

#                 for trackname, track_hist_func in track_hist_funcs.iteritems():
#                     # col_index is the index of the trackname
#                     #   (if 5 tracks, it can be 0, 1, 2, 3 or 4)
#                     col_index = chromosome.index_continuous(trackname)

#                     # select track data within the range of the segments and
#                     #   for the specified track
#                     continuous_col = continuous[start_i:end_i, col_index]

#                     # Mask continuous_col and label_map to only indices
#                     #   where the data is defined
#                     keep_mask = isfinite(continuous_col)
#                     continuous_col_ok = continuous_col[keep_mask]
#                     label_map_ok = label_map[keep_mask]

#                     # Keep track of number of defined data points
#                     if trackname in segmentation.tracks:
#                         nseg_dps += continuous_col_ok.shape[0]

#                     # This calls numpy.histogram2d
#                     hist, xedges, yedges = track_hist_func(label_map_ok,
#                                                            continuous_col_ok)
#                     res_trackname = res[trackname]
#                     assert (res_trackname[1] == xedges).all()
#                     assert (res_trackname[2] == yedges).all()
#                     res_trackname[0] += hist.astype(int)

#                     progress.next()

#             if quick: break

#     progress.end()
#     print ""
#     print "Segmentation tracks: %s" % str(segmentation.tracks)
#     print "Read %s non-NaN values from segmentation tracks\n" % nseg_dps
#     return res, nseg_dps

# def calc_stats2(histogram):
#     """Calculate approximate track statistics (mean, sd) for each label

#     histogram: the histogram returned by calc_histogram

#     Values are approximated from binned distributions
#     - mean is a lower-bound on the actual mean

#     Returns a dict: label_key -> dict( trackname -> {"mean", "sd", ...} )

#     """
#     stats = defaultdict(partial(defaultdict, dict))
#     for trackname, track_hist_data in histogram.iteritems():
#         (track_hist, label_edges, value_edges) = track_hist_data
#         for label_edge, hist_row in zip(label_edges, track_hist):
#             label_key = int(label_edge)
#             cur_stat = stats[trackname][label_key]
#             n = nansum(hist_row)
#             if n == 0:
#                 mean = NAN
#                 sd = NAN
#             else:
#                 mean = (hist_row * value_edges[:-1]).sum() / n
#                 sd = (hist_row * (value_edges[:1] - mean)**2).sum() / (n - 1)

#             cur_stat["sd"] = sd
#             cur_stat["mean"] = mean

#     return stats

# ## Saves the histogram data to a tab file
# def save_tab(labels, histogram, dirpath, clobber=False, group_labels=False,
#              basename=NAMEBASE, fieldnames=FIELDNAMES):
#     with tab_saver(dirpath, basename, fieldnames, clobber=clobber) as saver:
#         for label_key, label_histogram in histogram.iteritems():
#             label = labels[label_key]
#             for trackname, (hist, edges) in label_histogram.iteritems():
#                 for lower_edge, count in zip(edges, hist.tolist() + ["NA"]):
#                     saver.writerow(locals())
#             if group_labels: return  # Only write grouped label


# ## Saves the track stats to a tab file
# def save_stats_tab2(labels, stats, dirpath, clobber=False, group_labels=False,
#                    basename=NAMEBASE_STATS, fieldnames=FIELDNAMES_STATS):
#     with tab_saver(dirpath, basename, fieldnames, clobber=clobber) as saver:
#         for trackname, track_stats in stats.iteritems():
#             for label_key, label_stats in track_stats.iteritems():
#                 label = labels[label_key]
#                 mean = label_stats["mean"]
#                 sd = label_stats["sd"]
#                 saver.writerow(locals())

## Plots the tab file data and saves the plots
def save_plot(ntracks, nlabels, segtracks, dirpath,
              clobber=False, group_labels=False,
              basename=NAMEBASE, ecdf=False, mnemonics=[]):
    """
    ntracks: number of tracks in the Genome
    nlabels: number of labels in the Segmentation
    segtracks: a list of strings; the name of each track in the segmentation
    """
    start_R()

    if group_labels:  # Only one plot per track
        nlabels = 0
        mnemonics = []  # Ignore mnemonics

    # load data from the file just written out
    # it's bad to remake the filename
    # XXX: rewriting everything in terms of a class would allow avoiding this
    tabfilename = make_tabfilename(dirpath, basename)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    print >>sys.stderr, "\tPlotting distributions from: %s" % tabfilename
    if ecdf:
        mode = "ecdf"
    else:
        mode = "histogram"

    # XXX: Figure out ntracks and nlabels from file
    with image_saver(dirpath, basename,
                     width=PNG_BASE_WIDTH + \
                         ntracks * PNG_WIDTH_PER_SUBPLOT,
                     height=PNG_BASE_HEIGHT + \
                         (nlabels + 1) * PNG_HEIGHT_PER_SUBPLOT,
                     clobber=clobber,
                     downsample=True):  # Speeds up plotting
        r.plot(r["plot.signal"](tabfilename, segtracks=segtracks,
                                mnemonics=mnemonics, mode=mode))

def save_stats_plot(dirpath, basename=NAMEBASE_STATS, datafilename=None,
                    clobber=False, mnemonics=[], translation_file=None,
                    gmtk=False):
    """
    if datafilename is specified, it overrides dirpath/basename.tab as
    the data file for plotting.

    """
    ## Plot the track stats data with R
    start_R()
    r_source("track_statistics.R")

    if datafilename is None:
        datafilename = make_tabfilename(dirpath, basename)

    if not os.path.isfile(datafilename):
        die("Unable to find stats data file: %s" % datafilename)

    if translation_file is None:
        translation_file = ""  # None is not currently supported in rpy2

    with image_saver(dirpath, basename,
                     width=STATS_PNG_SIZE,
                     height=STATS_PNG_SIZE,
                     clobber=clobber):
        r.plot(r["plot.track.stats"](datafilename, mnemonics=mnemonics,
                                     translation_file=translation_file,
                                     gmtk=gmtk))

def save_html(dirpath, genomedatadir, nseg_dps=None,
              ecdf=False, clobber=False):
    if ecdf:
        title = "%s (ECDF mode)" % HTML_TITLE
    else:
        title = HTML_TITLE

    title = "%s (%s)" % (title, os.path.basename(genomedatadir))

    if nseg_dps is None:
        nseg_dps = "???"

    extra_namebases = {"stats": NAMEBASE_STATS}
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, title=title, segdatapoints=nseg_dps,
                  extra_namebases=extra_namebases)

## Package entry point
def validate(bedfilename, genomedatadir, dirpath, group_labels=False,
             clobber=False, calc_ranges=False, inputdirs=None,
             quick=False, replot=False, noplot=False,
             nbins=NBINS, value_range=(None, None),
             ecdf=False, mnemonicfilename=None,
             create_mnemonics=False, chroms=None):

    setup_directory(dirpath)
    genome = Genome(genomedatadir)
    segmentation = Segmentation(bedfilename)

    ntracks = genome.num_tracks_continuous  # All tracks; not just segtracks
    segtracks = segmentation.tracks
    labels = segmentation.labels

    if group_labels:
        nbins = None  # Set bins automatically
        nlabels = 0  # Not using labels
        fieldnames = FIELDNAMES[1:]  # No 'label' field
    else:
        nlabels = len(labels)
        fieldnames = FIELDNAMES

    if inputdirs:
        histogram = SignalHistogram()
        for inputdir in inputdirs:
            try:
                sub_histogram = SignalHistogram.read(inputdir)
            except IOError, e:
                print >>sys.stderr, "Problem reading data from %s: %s" % \
                    (inputdir, e)
            else:
                histogram.add(sub_histogram)

    elif not replot:
        # Generate histogram
        histogram = SignalHistogram.calculate(genome, segmentation,
                                              nbins=nbins,
                                              group_labels=group_labels,
                                              calc_ranges=calc_ranges,
                                              value_range=value_range,
                                              quick=quick, chroms=chroms)
    else:
        histogram = None

    if not replot:
        histogram.save(dirpath, clobber=clobber)

        stats = calc_stats(histogram)
        save_stats_tab(stats, dirpath, clobber=clobber,
                       basename=NAMEBASE_STATS)

        if mnemonicfilename is None and create_mnemonics:
            statsfilename = make_tabfilename(dirpath, NAMEBASE_STATS)
            mnemonicfilename = create_mnemonic_file(statsfilename, dirpath,
                                                    clobber=clobber)

    mnemonics = map_mnemonics(segmentation.labels, mnemonicfilename)

    if not noplot:
        save_stats_plot(dirpath, basename=NAMEBASE_STATS,
                        clobber=clobber, mnemonics=mnemonics)
        #save_plot(ntracks, nlabels, segtracks, dirpath, clobber=clobber,
        #          group_labels=group_labels, ecdf=ecdf, mnemonics=mnemonics)

    if histogram:
        try:
            nseg_dps = histogram.metadata["nseg_dps"]
        except KeyError:
            nseg_dps = None
        save_html(dirpath, genomedatadir, nseg_dps=nseg_dps,
                  ecdf=ecdf, clobber=clobber)


def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] BEDFILE GENOMEDATADIR"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    group.add_option("--clobber", action="store_true",
                     dest="clobber", default=False,
                     help="Overwrite existing output files if the specified"
                     " directory already exists.")
    group.add_option("--quick", action="store_true",
                     dest="quick", default=False,
                     help="Compute values only for one chromosome.")
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate plots")
    group.add_option("--replot", action="store_true",
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--create-mnemonics", action="store_true",
                     dest="create_mnemonics", default=False,
                     help="If mnemonics are not specified, they will be"
                     " created and used for plotting")
#     group.add_option("--group-labels", action="store_true",
#                      dest="group_labels", default=False,
#                      help="Group track distributions over all labels."
#                      " BEDFILE will be ignored")
    group.add_option("--ecdf", action="store_true",
                     dest="ecdf", default=False,
                     help="Plot empiracle cumulative density inside each panel"
                     " instead of a normal histogram (turns off log-y)")
    group.add_option("--calc-ranges", action="store_true",
                     dest="calc_ranges", default=False,
                     help="Calculate ranges for distribution plots from"
                     " segmentation data (slower) instead of using whole"
                     " genome data (default).")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Histogram options")
    group.add_option("-n", "--nbins", type="int", metavar="N",
                     dest="nbins", default=NBINS,
                     help="Number of bins for signal distribution"
                     " [default: %default]")
    group.add_option("-c", "--chrom", action="append", metavar="CHROM",
                     dest="chroms", default=None,
                     help="Only perform the analysis on data in CHROM,"
                     " where CHROM is a chromosome name such as chr21 or"
                     " chrX (option can be used multiple times to allow"
                     " multiple chromosomes)")
#     group.add_option("--min-value", type="float", metavar="N",
#                      dest="min_value", default=None,
#                      help="Minimum signal track value used in binning"
#                      " (overrides min from --calc-ranges)"
#                      " (values below will be ignored)")
#     group.add_option("--max-value", type="float", metavar="N",
#                      dest="max_value", default=None,
#                      help="Maximum signal track value used in binning"
#                      " (overrides max from --calc-ranges)"
#                      " (values above will be ignored)")
    parser.add_option_group(group)

    group = OptionGroup(parser, "I/O options")
    group.add_option("--mnemonic-file", dest="mnemonicfilename",
                     default=None, metavar="FILE",
                     help="If specified, labels will be shown using"
                     " mnemonics found in FILE")
    group.add_option("-i", "--indir", dest="inputdirs",
                     default=None, action="append", metavar="DIR",
                     help="Load data from this directory"
                     " This directory should be the output"
                     " directory of a previous run of this module."
                     " This option can be specified multiple times to"
                     " merge previous results together.")
    group.add_option("-o", "--outdir", metavar="DIR",
                     dest="outdir", default="%s" % MODULE,
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    if options.inputdirs:
        for inputdir in options.inputdirs:
            if inputdir == options.outdir:
                parser.error("Output directory cannot be an input directory")

    if options.noplot and options.replot:
        parser.error("noplot and replot are contradictory")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]
    genomedatadir = args[1]
    kwargs = {"clobber": options.clobber,
#               "group_labels": options.group_labels,
              "quick": options.quick,
              "calc_ranges": options.calc_ranges,
              "replot": options.replot,
              "noplot": options.noplot,
              "nbins": options.nbins,
#               "value_range": (options.min_value, options.max_value),
              "ecdf": options.ecdf,
              "inputdirs": options.inputdirs,
              "mnemonicfilename": options.mnemonicfilename,
              "create_mnemonics": options.create_mnemonics,
              "chroms": options.chroms}

    validate(bedfilename, genomedatadir, options.outdir, **kwargs )

if __name__ == "__main__":
    sys.exit(main())
