#!/usr/bin/env python
from __future__ import division, with_statement

"""
Provides command-line and package entry points for analyzing the
length distribution of entries in a provided segmentation or annotation file.

"""

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE="length_distribution"

import os
import sys

from collections import defaultdict
from numpy import concatenate, median

from . import Annotation, die, RInterface, open_transcript, \
     add_common_options
from .common import get_ordered_labels, LABEL_ALL, make_tabfilename, \
     setup_directory, tab_saver

from .html import save_html_div
from .version import __version__

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

R = RInterface(["common.R", "length.R"])

## Given an annotation file, returns the length of each entry
def annotation_lengths(annotation):
    # key: label_key
    # val: list of numpy.ndarray
    lengths = defaultdict(list)
    labels = annotation.labels

    # convert segment coords to lengths
    for rows in annotation.chromosomes.itervalues():
        for label_key in labels.iterkeys():
            labeled_row_indexes = (rows['key'] == label_key)
            labeled_rows = rows[labeled_row_indexes]
            length = labeled_rows['end'] - labeled_rows['start']
            lengths[label_key].append(length)
            lengths[LABEL_ALL].append(length)

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
        total_bp = num_bp[LABEL_ALL]
        row = make_size_row(LABEL_ALL, lengths[LABEL_ALL],
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
              mnemonic_file=None, transcriptfile=None, ropts=None):
    R.start(transcriptfile=transcriptfile, verbose=verbose)

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.length", dirpath, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber, ropts=ropts)

## Generates and saves an R plot of the length distributions
def save_size_plot(dirpath, namebase=NAMEBASE_SIZES, clobber=False,
                   verbose=True, mnemonic_file=None,
                   show_bases=True, show_segments=True,
                   transcriptfile=None, ropts=None):
    R.start(transcriptfile=transcriptfile, verbose=verbose)

    # Load data from corresponding tab file
    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.segment.sizes", dirpath, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber,
           show_segments=show_segments, show_bases=show_bases, ropts=ropts)

def save_html(dirpath, clobber=False, mnemonicfile=None, verbose=True):
    extra_namebases = {"sizes": NAMEBASE_SIZES}
    save_html_div(TEMPLATE_FILENAME, dirpath, namebase=NAMEBASE,
                  extra_namebases=extra_namebases, mnemonicfile=mnemonicfile,
                  tables={"": NAMEBASE_SIZES}, module=MODULE,
                  clobber=clobber, title=HTML_TITLE, verbose=verbose)

## Package entry point
def validate(filename, dirpath, clobber=False, replot=False, noplot=False,
             verbose=True, mnemonic_file=None, ropts=None,
             show_segments=True, show_bases=True):
    if not replot:
        setup_directory(dirpath)
        annotation = Annotation(filename, verbose=verbose)

        labels = annotation.labels

        lengths, num_bp = annotation_lengths(annotation)
        save_tab(lengths, labels, num_bp, dirpath,
                 clobber=clobber, verbose=verbose)
        save_size_tab(lengths, labels, num_bp, dirpath,
                      clobber=clobber, verbose=verbose)

    if not noplot:
        with open_transcript(dirpath, MODULE) as transcriptfile:
            save_plot(dirpath, mnemonic_file=mnemonic_file,
                      clobber=clobber, verbose=verbose,
                      transcriptfile=transcriptfile, ropts=ropts)
            save_size_plot(dirpath, clobber=clobber, verbose=verbose,
                           mnemonic_file=mnemonic_file,
                           show_segments=show_segments,
                           show_bases=show_bases, ropts=ropts,
                           transcriptfile=transcriptfile)

    save_html(dirpath, clobber=clobber, mnemonicfile=mnemonic_file,
              verbose=verbose)

def parse_options(args):
    from optparse import OptionGroup, OptionParser

    usage = "%prog [OPTIONS] ANNOTATION"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'replot', 'noplot'])
    group.add_option("--no-segments", action="store_false",
                     dest="show_segments", default=True,
                     help="Do not show total segments covered by labels"
                     " in size plot")
    group.add_option("--no-bases", action="store_false",
                     dest="show_bases", default=True,
                     help="Do not show total bases covered by labels"
                     " in size plot")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    filename = args[0]

    kwargs = dict(options.__dict__)
    outdir = kwargs.pop("outdir")
    validate(filename, outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
