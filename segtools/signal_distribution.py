#!/usr/bin/env python
from __future__ import division, with_statement

"""
Provides command-line and package entry points for analyzing the
distribution of signal values conditioned on segment labels.
"""

# A package-unique, descriptive module name used for filenames, etc
MODULE="signal_distribution"


import os
import sys

from collections import defaultdict
from genomedata import Genome
from functools import partial
from numpy import arcsinh, isfinite, isnan, longdouble, sqrt, square, zeros

from . import log, Segmentation, die, RInterface, add_common_options, \
     open_transcript, ProgressBar
from .common import iter_segments_continuous,  \
     make_tabfilename, setup_directory, tab_reader, tab_saver
from .html import save_html_div
from .mnemonics import create_mnemonic_file
from .version import __version__

FIELDNAMES = ["label", "trackname", "mean", "sd", "n"]
NAMEBASE = str(MODULE)

HTML_TITLE = "Signal value statistics"
HTML_TEMPLATE_FILENAME = "signal_div.tmpl"

R = RInterface(["common.R", "signal.R"])

class SignalStats(object):
    def __init__(self, data=None, **fields):
        self._data = data
        self._metadata = fields

    @property
    def data(self):
        return self._data

    @property
    def metadata(self):
        return self._metadata

    @staticmethod
    def from_segmentation(genome, segmentation, chroms=None, transformation=None,
                          quick=False, verbose=True):
        """
        Computes a mean and variance for each track-label pair.

        if chroms is:
        - a sequence of chromosome names: only those chromosomes are processed
        - False, [], None: all chromosomes are processed

        if quick is:
        - True: stats are calculated for only the first chromosome

        returns SignalStats object from data
        """
        assert genome is not None
        assert segmentation is not None

        labels = segmentation.labels
        tracks = genome.tracknames_continuous
        # A dict from tracks to a range tuple
        track_indices = dict(zip(tracks, range(len(tracks))))

        if len(tracks) == 0:
            die("Trying to calculate histogram for no tracks")

        (sum_total, sum2_total, dp_total) = \
            (zeros((len(tracks), len(labels)), dtype=longdouble),
             zeros((len(tracks), len(labels)), dtype=longdouble),
             zeros((len(tracks), len(labels)), dtype=int))
        log("Generating signal distribution histograms", verbose)

        with genome:
            if chroms:
                chromosomes = [genome[chrom] for chrom in chroms]
            else:
                chromosomes = genome

            for chromosome in chromosomes:
                chrom = chromosome.name
                if verbose:
                    try:
                        segments = segmentation.chromosomes[chrom]
                    except KeyError:
                        continue
                    progress = ProgressBar(len(segments) * len(tracks),
                                           label="  %s: " % chrom)

                for track_i, track in enumerate(tracks):
                    col_index = chromosome.index_continuous(track)
                    col_sum = sum_total[col_index]
                    col_sum2 = sum2_total[col_index]
                    col_dp = dp_total[col_index]
                    # Iterate through supercontigs and segments together
                    for segment, seg_data in \
                            iter_segments_continuous(chromosome, segmentation,
                                                     column=col_index,
                                                     verbose=verbose):
                        seg_label = segment['key']
                        seg_data_nonan = seg_data[isfinite(seg_data)]
                        if transformation == "arcsinh":
                            seg_data_nonan = arcsinh((seg_data_nonan))


                        col_sum[seg_label] += seg_data_nonan.sum(dtype=longdouble)
                        col_sum2[seg_label] += square(seg_data_nonan).sum(dtype=longdouble)
                        col_dp[seg_label] += seg_data_nonan.shape[0]
                        if verbose:
                            progress.next()

                progress.end()
                if quick: break  # 1 chromosome

        means = sum_total / dp_total
        sds = sqrt((sum2_total - (square(sum_total) / dp_total))/(dp_total - 1))

        stats = defaultdict(partial(defaultdict, dict))
        for label in labels:
            for trackname, track_index in track_indices.iteritems():
                # cur_stat = stats[label][trackname]
                # Paul:
                # Use the name of the label name defined in the
                # segmentation file as the key to the stats dictionary
                # so segtools internal labeling is not propagated to
                # user output.
                segmentfile_defined_label = labels[label]
                cur_stat = stats[segmentfile_defined_label][trackname]
                cur_stat["sd"] = sds[track_index, label]
                cur_stat["mean"] = means[track_index, label]
                cur_stat["n"] = dp_total[track_index, label]

        return SignalStats(stats)

    @staticmethod
    def from_file(dirpath, namebase=NAMEBASE, fieldnames=FIELDNAMES, **kwargs):
        """Load statistics from the tab file in an output directory"""
        stats = defaultdict(partial(defaultdict, dict))
        with tab_reader(dirpath, namebase, fieldnames=fieldnames,
                        **kwargs) as (reader, metadata):
            for row in reader:
                label = row.pop('label')
                trackname = row.pop('trackname')
                stats[label][trackname] = row

        return SignalStats(stats, **metadata)

    def add(self, o_stats):
        """Add data from a second SignalStats object into this one"""
        assert isinstance(o_stats, SignalStats)
        if self.data is None:
            self._data = o_stats.data
        else:
            for label, label_stats in self._data.iteritems():
                for trackname, track_stats in label_stats.iteritems():
                    o_track_stats = o_stats.data[label][trackname]
                    n_1 = int(track_stats['n'])
                    n_2 = int(o_track_stats['n'])
                    n = n_1 + n_2

                    s_1 = float(track_stats['mean']) * n_1
                    s_2 = float(o_track_stats['mean']) * n_2
                    s = s_1 + s_2

                    s2_1 = (s_1 - float(track_stats['sd']) * (n_1 - 1)) * n_1
                    s2_2 = (s_2 - float(o_track_stats['sd']) * (n_2 - 1)) * n_2
                    s2 = s2_1 + s2_2

                    track_stats['mean'] = s / n
                    track_stats['sd'] = (s2 - (square(s) / n))/(n - 1)
                    track_stats['n'] = n

    def save_tab(self, dirpath, clobber=False, verbose=True,
                 namebase=NAMEBASE, fieldnames=FIELDNAMES):
        with tab_saver(dirpath, namebase, fieldnames, verbose=verbose,
                       clobber=clobber) as saver:
            for label, label_stats in self._data.iteritems():
                for trackname, track_stats in label_stats.iteritems():
                    # assignments used for locals()
                    mean = track_stats["mean"]
                    sd = track_stats["sd"]
                    n = track_stats["n"]
                    saver.writerow(locals())

    @staticmethod
    def save_plot(dirpath, namebase=NAMEBASE, filename=None,
                  clobber=False, mnemonic_file=None, translation_file=None,
                  allow_regex=False, gmtk=False, verbose=True,
                  label_order=[], track_order=[], ropts=None,
                  transcriptfile=None):
        """
        if filename is specified, it overrides dirpath/namebase.tab as
        the data file for plotting.
        """
        ## Plot the track stats data with R
        R.start(verbose=verbose, transcriptfile=transcriptfile)
        R.source("track_statistics.R")

        if filename is None:
            filename = make_tabfilename(dirpath, namebase)

        if not os.path.isfile(filename):
            die("Unable to find stats data file: %s" % filename)

        R.plot("save.track.stats", dirpath, namebase, filename,
               mnemonic_file=mnemonic_file,
               translation_file=translation_file, ropts=ropts,
               as_regex=allow_regex, gmtk=gmtk, clobber=clobber,
               label_order=label_order, track_order=track_order)


def save_html(dirpath, genomedatadir, verbose=True, clobber=False):
    title = "%s (%s)" % (HTML_TITLE, os.path.basename(genomedatadir))

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, title=title, verbose=verbose)

def read_order_file(filename):
    if filename is None:
        return []
    elif os.path.isfile(filename):
        order = []
        with open(filename) as ifp:
            for line in ifp:
                line = line.strip()
                if line:
                    order.append(line)

        return order
    else:
        raise IOError("Could not find order file: %s" % filename)


## Package entry point
def validate(bedfilename, genomedatadir, dirpath, clobber=False,
             quick=False, replot=False, noplot=False, verbose=True,
             mnemonic_file=None, create_mnemonics=False,
             inputdirs=None, chroms=None, ropts=None,
             label_order_file=None, track_order_file=None, transformation=None):

    if not replot:
        setup_directory(dirpath)
        genome = Genome(genomedatadir)
        segmentation = Segmentation(bedfilename, verbose=verbose)

    if inputdirs:
        # Merge stats from many input directories
        stats = SignalStats()
        for inputdir in inputdirs:
            try:
                sub_stats = SignalStats.from_file(inputdir, verbose=verbose)
            except IOError, e:
                log("Problem reading data from %s: %s" % (inputdir, e))
            else:
                stats.add(sub_stats)
    elif replot:
        stats = SignalStats.from_file(dirpath, verbose=verbose)
    else:
        # Calculate stats over segmentation
        stats = SignalStats.from_segmentation(genome, segmentation, transformation=transformation,
                                              quick=quick, chroms=chroms,
                                              verbose=verbose)

    if not replot:
        stats.save_tab(dirpath, clobber=clobber, verbose=verbose)

        if mnemonic_file is None and create_mnemonics:
            statsfilename = make_tabfilename(dirpath, NAMEBASE)
            mnemonic_file = create_mnemonic_file(statsfilename, dirpath,
                                                 clobber=clobber,
                                                 verbose=verbose)

    if not noplot:
        if label_order_file is not None:
            log("Reading label ordering from: %s" % label_order_file)
        label_order = read_order_file(label_order_file)

        if track_order_file is not None:
            log("Reading track ordering from: %s" % track_order_file)
        track_order = read_order_file(track_order_file)

        with open_transcript(dirpath, MODULE) as transcriptfile:
            stats.save_plot(dirpath, namebase=NAMEBASE, clobber=clobber,
                            mnemonic_file=mnemonic_file, verbose=verbose,
                            label_order=label_order, track_order=track_order,
                            ropts=ropts, transcriptfile=transcriptfile)

    save_html(dirpath, genomedatadir, clobber=clobber, verbose=verbose)


def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] SEGMENTATION GENOMEDATAFILE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'quick', 'noplot',
                               'replot'])
    group.add_option("--create-mnemonics", action="store_true",
                     dest="create_mnemonics", default=False,
                     help="If mnemonics are not specified, they will be"
                     " created and used for plotting")
    parser.add_option_group(group)

    group = OptionGroup(parser, "I/O options")
    group.add_option("-t", "--transformation", action="store", metavar="TRANSFORMATION",
                     dest="transformation", default=None,
                     help="Applies the transformation on the data upon reading  from genomedata."
					 "The default is None, and currently only 'arcsinh' is implemented.")
    group.add_option("-c", "--chrom", action="append", metavar="CHROM",
                     dest="chroms", default=None,
                     help="Only perform the analysis on data in CHROM,"
                     " where CHROM is a chromosome name such as chr21 or"
                     " chrX (option can be used multiple times to allow"
                     " multiple chromosomes)")
    group.add_option("--order-tracks", dest="track_order_file",
                     default=None, metavar="FILE",
                     help="If specified, tracks will be displayed"
                     " in the order in FILE. FILE must be a permutation"
                     " of all the printed tracks, one per line, exact"
                     " matches only.")
    group.add_option("--order-labels", dest="label_order_file",
                     default=None, metavar="FILE",
                     help="If specified, label will be displayed"
                     " in the order in FILE. FILE must be a permutation"
                     " of all the labels (after substituting with mnemonics,"
                     " if specified), one per line, exact"
                     " matches only.")
    group.add_option("-i", "--indir", dest="inputdirs",
                     default=None, action="append", metavar="DIR",
                     help="Load data from this directory."
                     " This directory should be the output"
                     " directory of a previous run of this module."
                     " This option can be specified multiple times to"
                     " merge previous results together.")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    if options.inputdirs:
        for inputdir in options.inputdirs:
            if inputdir == options.outdir:
                parser.error("Output directory cannot be an input directory")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]
    genomedatadir = args[1]

    kwargs = dict(options.__dict__)
    outdir = kwargs.pop('outdir')
    validate(bedfilename, genomedatadir, outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
