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
#from contextlib import contextmanager
from itertools import repeat
from functools import partial
from numpy import array, concatenate, histogram, isfinite, nanmax, nanmin, \
    NINF, PINF, where
from rpy2.robjects import r, numpy2ri

# XXX: Do this without the kludgy constants
from .common import die, fill_array, \
    get_supercontig_splice, \
    image_saver, LABEL_ALL, load_segmentation, load_genome, \
    loop_segments_continuous, make_tabfilename, \
    map_mnemonics, map_segments, r_source, \
    SEGMENT_START_COL, SEGMENT_END_COL, SEGMENT_LABEL_KEY_COL, \
    setup_directory, tab_saver

from .html import save_html_div

FIELDNAMES = ["label", "trackname", "lower_edge", "count"]
NAMEBASE = "%s" % MODULE

HTML_TITLE = "Signal value distribution"
HTML_TEMPLATE_FILENAME = "signal_div.tmpl"

PNG_WIDTH_PER_SUBPLOT = 200
PNG_HEIGHT_PER_SUBPLOT = 150
PNG_BASE_HEIGHT = 200
PNG_BASE_WIDTH = 100

NUM_BINS = 100


def constant_factory(val):
    return repeat(val).next

def remove_nans(numbers):
# by Mirela
    return numbers[isfinite(numbers)]

## Returns a dict of trackname -> tuples: (min, max), one for each trackname
def seg_min_max(chromosome, segmentation, verbose=False):
    #print >>sys.stderr, "Looking for min & max in chromosome ", chromosome.name, " for track ", track_index
    #print >>sys.stderr, "segmentation: ", segmentation.chromosomes
    
    tracknames = chromosome.tracknames_continuous
    limits = dict([(trackname, (PINF, NINF)) for trackname in tracknames])
    
    if chromosome.name not in segmentation.chromosomes.keys():
        return limits
    
    for segment, continuous in \
            loop_segments_continuous(chromosome, segmentation,
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

def get_num_tracks(genome):    
    # MA: this should really be part of the Genome class - TODO
    for chromosome in genome:    
        tracknames = chromosome.tracknames_continuous
        return len(tracknames)

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
                #print >>sys.stderr, "Limits for %s: %s" % (trackname, res[trackname]) 
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

# def print_array(tag, arr, type="%d"):
#     fstring = "%%s:  %s, %s, ..., %s, ..., %s, %s" % tuple([type]*5)
#     print fstring % (tag,
#                      arr[0],
#                      arr[1],
#                      arr[int(arr.shape[0] / 2)],
#                      arr[-2],
#                      arr[-1])

## Computes a set of histograms, one for each track-label pair.
def calc_histogram(genome, segmentation, num_bins=None,
                   group_labels=False, calc_ranges=False,
                   value_range=(None, None), quick=False):
    """
    if num_bins is:
    - None: There will be a bin for every possible value
    - > 0: There will be this many bins
      
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
      key: label_key
      val: dict
           key: trackname
           val: [hist, edges]
    """
    assert genome is not None
    assert segmentation is not None

    if calc_ranges and any(lim is None for lim in value_range):
        # Calculate ranges over segmentation only
        print >>sys.stderr, "Calculating track ranges...",
        track_ranges = load_track_ranges(genome, segmentation)
        print >>sys.stderr, "done."
    else:
        # Lookup ranges in genomedata
        track_ranges = load_track_ranges(genome)
        
    # Override ranges as necessary
    value_min, value_max = value_range
    for trackname, track_range in track_ranges.iteritems():
        track_min, track_max = track_range
        if value_min is not None:
            track_min = value_min
        if value_max is not None:
            track_max = value_max
        track_ranges[trackname] = (track_min, track_max)

    tracks = track_ranges.keys()  # Re-select just the tracks used

    if len(tracks) == 0:
        die("Trying to calculate histogram for no tracks")

        
    # Calculate bins for each track
    if num_bins is None:  # Have a bin for every possible integer value
        track_bins = dict([(trackname,
                            math.ceil(track_range[1]) -
                            math.floor(track_range[0]))
                           for trackname, track_range in track_ranges])
    else:  # Use the given number of bins
        assert num_bins > 0
        track_bins = dict([(trackname, num_bins) for trackname in tracks])

        
    # key: trackname
    # val: histogram_func: a function that generates a histogram from data
    #                      with uniform bins within a track
    # For every trackname, there's a function that generates histograms
    histogram_funcs = dict((trackname, partial(histogram,
                                               bins=track_bins[trackname],
                                               range=track_ranges[trackname],
                                               new=True))
                           for trackname in tracks)
    # Trackname is one the track names, track_range is (min,max)
    # functools.partial(func[, *args][, **keywords])
    #   Returns a new partial object which when called will behave like func
    #   called with args and keywords.
    
    labels = segmentation.labels
    if group_labels:  # Pretend there is just one "all" label
        labels = dict([(-1, "all")])
        
    # Return dictionary: see function docstring for specification
    # Now for every label and for every trackname
    #   there's a histogram array of the right boundary and frequency
    res = dict((label_key, dict((trackname,
                                 list(histogram_func(array([]))))
                                for trackname, histogram_func
                                in histogram_funcs.iteritems()))
               for label_key in labels)
    
    print >>sys.stderr, "Generating signal distribution histograms"
    num_tracks = len(tracks)
    num_labels = len(labels)
    for chromosome in genome:
        chrom = chromosome.name
        print >>sys.stderr, "\t%s" % chrom

        # Iterate through supercontigs and segments together
        for segment, continuous_seg in \
                loop_segments_continuous(chromosome, segmentation):

            if group_labels:
                cur_res = res  # All labels are used
            else:
                seg_label_key = segment[SEGMENT_LABEL_KEY_COL]
                cur_res = res[seg_label_key]
                
            # Iterate through each track
            for trackname, histogram_func in histogram_funcs.iteritems():
                # col_index is the index of the trackname
                #   (if 5 tracks, it can be 0, 1, 2, 3 or 4)
                col_index = chromosome.index_continuous(trackname)
                
                # continuous_col is the "intensity"
                #   (continuous number) from the data track
                # len(supercontig_map) = len(continuous_col)
                #   because they are for the same segment
                continuous_col = continuous_seg[:, col_index]
                
                # Remove the NaN's, otherwise the numpy.histogram function 
                #   gets confused and misses some of the data 
                cur_col_nonan = remove_nans(continuous_col)
                
                # Edges is a list containing the left edge of every bin,
                #   hist is the frequency for every bin
                # This calls numpy.histogram
                hist, edges = histogram_func(cur_col_nonan)
                
                res_trackname = cur_res[trackname]
                
                assert (res_trackname[1] == edges).all()
                res_trackname[0] += hist
        if quick: break  # 1 chromosome
    return res

## Saves the histogram data to a tab file
def save_tab(labels, histogram, dirpath, clobber=False, group_labels=False,
             namebase=NAMEBASE, fieldnames=FIELDNAMES):
    with tab_saver(dirpath, namebase, fieldnames, clobber=clobber) as saver:
        for label_key, label_histogram in histogram.iteritems():
            label = labels[label_key]
            for trackname, (hist, edges) in label_histogram.iteritems():
                for lower_edge, count in zip(edges, hist.tolist() + ["NA"]):
                    saver.writerow(locals())
            if group_labels: return  # Only write grouped label

## Plots the tab file data and saves the plots
def save_plot(num_tracks, num_labels, segtracks, dirpath,
              clobber=False, group_labels=False,
              basename=NAMEBASE, ecdf=False, mnemonics=[]):
    """
    num_tracks: number of tracks in the Genome
    num_labels: number of labels in the Segmentation
    segtracks: a list of strings; the name of each track in the segmentation
    """
    r_source("common.R")
    r_source("signal.R")

    if group_labels:  # Only one plot per track
        num_labels = 0
        mnemonics = []  # Ignore mnemonics

    # load data from the file just written out
    # it's bad to remake the filename
    # XXX: rewriting everything in terms of a class would allow avoiding this
    tabfilename = make_tabfilename(dirpath, basename)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    print >>sys.stderr, "\tPlotting distributions from: %s" % tabfilename
    # XXX: Figure out num_tracks and num_labels from file
    with image_saver(dirpath, basename,
                     width=PNG_BASE_WIDTH + \
                         num_tracks * PNG_WIDTH_PER_SUBPLOT,
                     height=PNG_BASE_HEIGHT + \
                         (num_labels + 1) * PNG_HEIGHT_PER_SUBPLOT,
                     clobber=clobber,
                     downsample=True):  # Speeds up plotting
        r.plot(r["plot.signal"](tabfilename, segtracks=segtracks,
                                mnemonics=mnemonics, ecdf=ecdf))

def save_html(dirpath, ecdf=False, clobber=False):
    if ecdf:
        title = "%s (ECDF mode)" % HTML_TITLE
    else:
        title = HTML_TITLE
        
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, title=title)
    
## Package entry point
def validate(bedfilename, genomedatadir, dirpath, group_labels=False,
             clobber=False, calc_ranges=False,
             quick=False, replot=False, noplot=False,
             num_bins=NUM_BINS, value_range=(None, None),
             ecdf=False, mnemonicfilename=None):
    
    setup_directory(dirpath)
    genome = load_genome(genomedatadir)
    segmentation = load_segmentation(bedfilename)
    
    assert genome is not None
    assert segmentation is not None
    
    mnemonics = map_mnemonics(segmentation.labels, mnemonicfilename)
    num_tracks = get_num_tracks(genome)  # All tracks (not just segtracks)
    segtracks = segmentation.tracks
    labels = segmentation.labels
    
    if group_labels:
        num_bins = None  # Set bins automatically
        num_labels = 0  # Not using labels
        fieldnames = FIELDNAMES[1:]  # No 'label' field
    else:
        num_labels = len(labels)
        fieldnames = FIELDNAMES
        
    if not replot:
        # Generate histogram and save tab and plot files
        histogram = calc_histogram(genome, segmentation, num_bins=num_bins,
                                   group_labels=group_labels,
                                   calc_ranges=calc_ranges,
                                   value_range=value_range,
                                   quick=quick)
    
        save_tab(labels, histogram, dirpath, clobber=clobber,
                 group_labels=group_labels, fieldnames=fieldnames)

    if not noplot:
        save_plot(num_tracks, num_labels, segtracks, dirpath, clobber=clobber,
                  group_labels=group_labels, ecdf=ecdf, mnemonics=mnemonics)

    save_html(dirpath, ecdf=ecdf, clobber=clobber)
                
            
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
    group.add_option("--replot", action="store_true", 
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--noplot", action="store_true", 
                     dest="noplot", default=False,
                     help="Do not generate plots")
    group.add_option("--group-labels", action="store_true", 
                     dest="group_labels", default=False,
                     help="Group track distributions over all labels."
                     " BEDFILE will be ignored")
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
    group.add_option("-n", "--num-bins", type="int",
                     dest="num_bins", default=NUM_BINS,
                     help="Number of bins for signal distribution"
                     " [default: %default]")
    group.add_option("--min-value", type="float",
                     dest="min_value", default=None,
                     help="Minimum signal track value used in binning"
                     " (overrides min from --calc-ranges)"
                     " (values below will be ignored)")
    group.add_option("--max-value", type="float",
                     dest="max_value", default=None,
                     help="Maximum signal track value used in binning"
                     " (overrides max from --calc-ranges)"
                     " (values above will be ignored)")
    parser.add_option_group(group)
    
    group = OptionGroup(parser, "I/O options")
    group.add_option("--mnemonic-file", dest="mnemonicfilename",
                     default=None,
                     help="If specified, labels will be shown using"
                     " mnemonics found in this file")
    group.add_option("-o", "--outdir", 
                     dest="outdir", default="%s" % MODULE,
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")
        
    if options.noplot and options.replot:
        parser.error("noplot and replot are contradictory")
        
    return (options, args)
    
## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]
    genomedatadir = args[1]
    value_range = (options.min_value, options.max_value)
    validate(bedfilename, genomedatadir, options.outdir, options.group_labels,
             clobber=options.clobber, quick=options.quick,
             calc_ranges=options.calc_ranges, replot=options.replot,
             noplot=options.noplot, num_bins=options.num_bins,
             value_range=value_range, ecdf=options.ecdf, 
             mnemonicfilename=options.mnemonicfilename)

if __name__ == "__main__":
    sys.exit(main())
