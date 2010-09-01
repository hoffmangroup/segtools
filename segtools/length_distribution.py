#!/usr/bin/env python
from __future__ import division, with_statement

"""
Provides command-line and package entry points for analyzing the segment
length distribution in a provided BED-formatted segmentation.

"""

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE="length_distribution"

__version__ = "$Revision$"

import os
import sys

from collections import defaultdict
from numpy import concatenate, median

from . import Segmentation
from .common import die, get_ordered_labels, LABEL_ALL, make_tabfilename, r_plot, r_source, setup_directory, tab_saver

from .html import save_html_div

FIELDNAMES_SUMMARY = ["label", "num.segs", "mean.len", "median.len",
                     "stdev.len", "num.bp", "frac.bp"]
FIELDNAMES = ["label", "length"]
NAMEBASE = "%s" % MODULE
NAMEBASE_SIZES = "segment_sizes"
TEMPLATE_FILENAME = "length_div.tmpl"

HTML_TITLE = "Length distribution"

PNG_WIDTH = 600
PNG_HEIGHT_BASE = 100  # Axes and label
PNG_HEIGHT_PER_LABEL = 45

def start_R():
    r_source("common.R")
    r_source("length.R")

## Given a segmentation, returns the length of each segment
def segmentation_lengths(segmentation):
    # key: label_key
    # val: list of numpy.ndarray
    lengths = defaultdict(list)
    labels = segmentation.labels

    # convert segment coords to lengths
    for segments in segmentation.chromosomes.itervalues():
        for label_key in labels.iterkeys():
            labeled_row_indexes = (segments['key'] == label_key)
            labeled_rows = segments[labeled_row_indexes]
            lengths[label_key].append(labeled_rows['end'] -
                                      labeled_rows['start'])

    # key: label_key
    # val: int
    num_bp = {}

    # convert lengths to:
    # key: label_key
    # val: numpy.ndarray(dtype=integer)
    for label_key, label_lengths_list in lengths.iteritems():
        label_lengths = concatenate(label_lengths_list)
        lengths[label_key] = label_lengths
        num_bp[label_key] = label_lengths.sum()

    return lengths, num_bp

## Specifies the format of a row of the summary tab file
def make_size_row(label, lengths, num_bp, total_bp):
    return {"label": label,
            "num.segs": len(lengths),
            "mean.len": "%.3f" % lengths.mean(),
            "median.len": "%.3f" % median(lengths),
            "stdev.len": "%.3f" % lengths.std(),
            "num.bp": num_bp,
            "frac.bp": "%.3f" % (num_bp / total_bp)}

## Saves the length summary data to a tab file, using mnemonics if specified
def save_size_tab(lengths, labels, num_bp, dirpath, verbose=True,
                  namebase=NAMEBASE_SIZES, clobber=False):
    ordered_keys, labels = get_ordered_labels(labels)
    with tab_saver(dirpath, namebase, FIELDNAMES_SUMMARY,
                   clobber=clobber, verbose=verbose) as saver:
        # "all" row first
        total_bp = sum(num_bp.itervalues())
        all_lengths = concatenate(lengths.values())
        row = make_size_row(LABEL_ALL, all_lengths,
                               total_bp, total_bp)
        saver.writerow(row)

        # Remaining label rows
        for label_key in ordered_keys:
            label = labels[label_key]
            label_lengths = lengths[label_key]
            num_bp_label = num_bp[label_key]
            row = make_size_row(label, label_lengths,
                                   num_bp_label, total_bp)
            saver.writerow(row)

def make_row(label, length):
    return {"label": label,
            "length": length}

def save_tab(lengths, labels, num_bp, dirpath, clobber=False, verbose=True):
    # fix order for iterating through dict; must be consistent
    label_keys = sorted(labels.keys())

    # XXXopt: allocating space in advance might be faster
    lengths_array = concatenate([lengths[label_key]
                                 for label_key in label_keys])

    # creates an array that has the label repeated as many times as
    # the length
    labels_array = concatenate([[labels[label_key]] * len(lengths[label_key])
                               for label_key in label_keys])

    with tab_saver(dirpath, NAMEBASE, FIELDNAMES, verbose=verbose,
                   clobber=clobber) as saver:
        for label, length in zip(labels_array, lengths_array):
            row = make_row(label, length)
            saver.writerow(row)

## Generates and saves an R plot of the length distributions
def save_plot(dirpath, namebase=NAMEBASE, clobber=False, verbose=True,
              mnemonic_file=None):
    start_R()

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    r_plot("save.length", dirpath, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber, verbose=verbose)

## Generates and saves an R plot of the length distributions
def save_size_plot(dirpath, namebase=NAMEBASE_SIZES, clobber=False,
                   verbose=True, mnemonic_file=None):
    start_R()

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    r_plot("save.segment.sizes", dirpath, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber, verbose=verbose)

def save_html(dirpath, clobber=False, mnemonicfile=None):
    extra_namebases = {"sizes": NAMEBASE_SIZES}
    save_html_div(TEMPLATE_FILENAME, dirpath, namebase=NAMEBASE,
                  extra_namebases=extra_namebases, mnemonicfile=mnemonicfile,
                  tables={"": NAMEBASE_SIZES}, module=MODULE,
                  clobber=clobber, title=HTML_TITLE)

## Package entry point
def validate(bedfilename, dirpath, clobber=False, replot=False, noplot=False,
             verbose=True, mnemonic_file=None):
    if not replot:
        setup_directory(dirpath)
        segmentation = Segmentation(bedfilename, verbose=verbose)

        labels = segmentation.labels

        lengths, num_bp=segmentation_lengths(segmentation)
        save_tab(lengths, labels, num_bp, dirpath,
                 clobber=clobber, verbose=verbose)
        save_size_tab(lengths, labels, num_bp, dirpath,
                      clobber=clobber, verbose=verbose)

    if not noplot:
        save_plot(dirpath, mnemonic_file=mnemonic_file,
                  clobber=clobber, verbose=verbose)
        save_size_plot(dirpath, clobber=clobber, verbose=verbose,
                       mnemonic_file=mnemonic_file)

    save_html(dirpath, clobber=clobber, mnemonicfile=mnemonic_file)

def parse_options(args):
    from optparse import OptionGroup, OptionParser

    usage = "%prog [OPTIONS] SEGMENTATION"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    group.add_option("--clobber", action="store_true",
                     dest="clobber", default=False,
                     help="Overwrite existing output files if the specified"
                     " directory already exists.")
    group.add_option("-q", "--quiet", action="store_false",
                     dest="verbose", default=True,
                     help="Do not print diagnostic messages.")
    group.add_option("--replot", action="store_true",
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate plots")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    group.add_option("-m", "--mnemonic-file", dest="mnemonic_file",
                     default=None, metavar="FILE",
                     help="If specified, labels will be shown using"
                     " mnemonics found in FILE.")
    group.add_option("-o", "--outdir",
                     dest="outdir", default="%s" % MODULE, metavar="DIR",
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]

    kwargs = {"clobber": options.clobber,
              "verbose": options.verbose,
              "replot": options.replot,
              "noplot": options.noplot,
              "mnemonic_file": options.mnemonic_file}
    validate(bedfilename, options.outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
