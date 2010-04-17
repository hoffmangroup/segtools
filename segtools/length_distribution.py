#!/usr/bin/env python
from __future__ import division, with_statement

"""
validate.py:

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
from rpy2.robjects import r, numpy2ri

from . import Segmentation
from .common import die, get_ordered_labels, LABEL_ALL, make_tabfilename, map_mnemonics, r_source, setup_directory, tab_saver

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
def save_size_tab(lengths, labels, num_bp, dirpath,
                  namebase=NAMEBASE_SIZES, clobber=False):
    ordered_keys, labels = get_ordered_labels(labels)
    with tab_saver(dirpath, namebase, FIELDNAMES_SUMMARY,
                   clobber=clobber) as saver:
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

def save_tab(lengths, labels, num_bp, dirpath, clobber=False):
    # fix order for iterating through dict; must be consistent
    label_keys = sorted(labels.keys())

    # XXXopt: allocating space in advance might be faster
    lengths_array = concatenate([lengths[label_key]
                                 for label_key in label_keys])

    # creates an array that has the label repeated as many times as
    # the length
    labels_array = concatenate([[labels[label_key]] * len(lengths[label_key])
                               for label_key in label_keys])

    with tab_saver(dirpath, NAMEBASE, FIELDNAMES, clobber=clobber) as saver:
        for label, length in zip(labels_array, lengths_array):
            row = make_row(label, length)
            saver.writerow(row)

## Generates and saves an R plot of the length distributions
def save_plot(dirpath, namebase=NAMEBASE, clobber=False,
              mnemonic_file=""):
    start_R()

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    if not mnemonic_file:
        mnemonic_file = ""  # None cannot be passed to R

    r["save.length"](dirpath, namebase, tabfilename,
                     mnemonic_file=mnemonic_file,
                     clobber=clobber)

## Generates and saves an R plot of the length distributions
def save_size_plot(dirpath, namebase=NAMEBASE_SIZES, clobber=False,
                   mnemonic_file=""):
    start_R()

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    if not mnemonic_file:
        mnemonic_file = ""  # None cannot be passed to R

    r["save.segment.sizes"](dirpath, namebase, tabfilename,
                            clobber=clobber)

def save_html(dirpath, clobber=False, mnemonicfile=None):
    extra_namebases = {"sizes": NAMEBASE_SIZES}
    save_html_div(TEMPLATE_FILENAME, dirpath, namebase=NAMEBASE,
                  extra_namebases=extra_namebases, mnemonicfile=mnemonicfile,
                  tables={"":NAMEBASE_SIZES}, module=MODULE,
                  clobber=clobber, title=HTML_TITLE)

## Package entry point
def validate(bedfilename, dirpath, clobber=False, replot=False, noplot=False,
             mnemonicfilename=None):
    if not replot:
        setup_directory(dirpath)
        segmentation = Segmentation(bedfilename)

        labels = segmentation.labels
        mnemonics = map_mnemonics(labels, mnemonicfilename)

        lengths, num_bp=segmentation_lengths(segmentation)
        save_tab(lengths, labels, num_bp, dirpath, clobber=clobber)
        save_size_tab(lengths, labels, num_bp, dirpath, clobber=clobber)

    if not noplot:
        save_plot(dirpath, clobber=clobber, mnemonic_file=mnemonicfilename)
        save_size_plot(dirpath, clobber=clobber,
                       mnemonic_file=mnemonicfilename)

    save_html(dirpath, clobber=clobber, mnemonicfile=mnemonicfilename)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] BEDFILE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("--clobber", action="store_true",
                      dest="clobber", default=False,
                      help="Overwrite existing output files if the specified"
                      " directory already exists.")
    parser.add_option("--replot", action="store_true",
                      dest="replot", default=False,
                      help="Load data from output tab files and"
                      " regenerate plots instead of recomputing data")
    parser.add_option("--noplot", action="store_true",
                      dest="noplot", default=False,
                      help="Do not generate plots")

    parser.add_option("--mnemonic-file", dest="mnemonicfilename",
                      default=None,
                      help="If specified, labels will be shown using"
                      " mnemonics found in this file")
    parser.add_option("-o", "--outdir",
                      dest="outdir", default="%s" % MODULE,
                      help="File output directory (will be created"
                      " if it does not exist) [default: %default]")

    (options, args) = parser.parse_args(args)

    if len(args) < 1:
        parser.error("Insufficient number of arguments")

    if options.noplot and options.replot:
        parser.error("noplot and replot are contradictory")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]

    validate(bedfilename, options.outdir, clobber=options.clobber,
             replot=options.replot, noplot=options.noplot,
             mnemonicfilename=options.mnemonicfilename)

if __name__ == "__main__":
    sys.exit(main())
