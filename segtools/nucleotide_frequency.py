#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
validate.py:

Provides command-line and package entry points for analyzing nucleotide
and dinucleotide frequencies for each segmentation label.

"""

import numpy
import os
import sys

from collections import defaultdict
from genomedata import Genome
from numpy import zeros
from rpy2.robjects import r, numpy2ri

from . import Segmentation
from .common import die, make_tabfilename, r_source, setup_directory, tab_saver

from .html import save_html_div

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE="nucleotide_frequency"

HTML_TITLE = "Nucleotide and dinucleotide content"
HTML_TEMPLATE_FILENAME = "nucleotide_div.tmpl"

# Fieldnames must agree with categories in order and content
# with a "label" category at the beginning and nuc before dinuc
FIELDNAMES = ["label", "A|T", "C|G", "?", "AA|TT", "AC|GT",
              "AG|CT", "AT", "CA|TG", "CC|GG", "CG", "GA|TC", "GC", "TA",
              "??"]
NAMEBASE = "%s" % MODULE

NUC_CATEGORIES = [('A','T'), ('C', 'G')]
DINUC_CATEGORIES = [[('A', 'A'), ('T', 'T')],
                    [('A', 'C'), ('G', 'T')],
                    [('A', 'G'), ('C', 'T')],
                    [('A', 'T')],
                    [('C', 'A'), ('T', 'G')],
                    [('C', 'C'), ('G', 'G')],
                    [('C', 'G')],
                    [('G', 'A'), ('T', 'C')],
                    [('G', 'C')],
                    [('T', 'A')]]

def start_R():
    r_source("common.R")
    r_source("dinucleotide.R")

## Caclulates nucleotide and dinucleotide frequencies over the specified
## segmentation
def calc_nucleotide_frequencies(segmentation, genome,
                                nuc_categories=NUC_CATEGORIES,
                                dinuc_categories=DINUC_CATEGORIES,
                                quick=False):

    # Store categories efficiently as dict
    # from each entry directly to category index
    quick_nuc_categories = defaultdict(dict)
    for index, category in enumerate(nuc_categories):
        for entry in category:
            #print("Mapping %s -> %d" % (entry, index))
            quick_nuc_categories[entry] = index

    quick_dinuc_categories = defaultdict(dict)
    for index, category in enumerate(dinuc_categories):
        for entry in category:
            #print("Mapping %s%s -> %d" % (entry[0], entry[1], index))
            quick_dinuc_categories[entry] = index

    # Store counts of each (di)nucleotide observed
    # separated by label
    nuc_counts = defaultdict(dict)
    dinuc_counts = defaultdict(dict)

    labels = segmentation.labels
    for label_key in labels.iterkeys():
        nuc_counts[label_key] = zeros(len(nuc_categories)+1, dtype=numpy.long)
        dinuc_counts[label_key] = zeros(len(dinuc_categories)+1,
                                        dtype=numpy.long)

    # Count (di)nucleotides over segmentation
    for chromosome in genome:
        chrom = chromosome.name
        print >>sys.stderr, "\t%s" % chrom

        try:
            segments = segmentation.chromosomes[chrom]
        except KeyError:
            continue

        for segment in segments:
            seg_start, seg_end, seg_label = segment

            sequence = chromosome.seq[seg_start:seg_end]
            sequence = sequence.tostring().upper()

            # XXXopt: could be optimized to use matrix operations
            # to determine number of occurances of each (di)nuc
            # Inch through iteration to look at pairs efficiently
            cur_nuc = None
            prev_nuc = None
            for nuc in sequence:
                prev_nuc = cur_nuc
                cur_nuc = nuc

                # Record nucleotide
                nuc_cat = quick_nuc_categories[cur_nuc]
                if nuc_cat == {}:  # Didn't find; put in last bin
                    #print >>sys.stderr, "Unknown nuc: %s" % (cur_nuc)
                    nuc_cat = len(nuc_counts[seg_label]) - 1

                nuc_counts[seg_label][nuc_cat] += 1

                # Record dinucleotide (last and current)
                if prev_nuc is not None:
                    dinuc_cat = quick_dinuc_categories[(prev_nuc, cur_nuc)]
                    if dinuc_cat == {}:  # Didn't find; put in last bin
                        dinuc_cat = len(dinuc_counts[seg_label]) - 1

                    dinuc_counts[seg_label][dinuc_cat] += 1

        # Only look at first chromosome if quick
        if quick: break

    return (nuc_counts, dinuc_counts)


def make_row(label, nuc_counts, dinuc_counts, fieldnames=FIELDNAMES):
    row = {}

    for i, field in enumerate(fieldnames):
        if field == "label":
            row[field] = label
        elif i <= len(nuc_counts):
            # Compensate for "label" field
            row[field] = nuc_counts[i - 1]
        else:
            # Compensate for "label" and nucleotide fields
            row[field] = dinuc_counts[i - len(nuc_counts) - 1]
    return row

# Fieldnames must agree with categories in order and content
# with a "label" category at the beginning, and nuc before dinuc
def save_tab(labels, nuc_counts, dinuc_counts, dirpath,
             fieldnames=FIELDNAMES, clobber=False):
    with tab_saver(dirpath, NAMEBASE, fieldnames, clobber=clobber) as saver:
        for label_key in sorted(labels.keys()):
            row = make_row(labels[label_key],
                           nuc_counts[label_key],
                           dinuc_counts[label_key])
            saver.writerow(row)

def save_plot(dirpath, clobber=False, mnemonic_file="", namebase=NAMEBASE):
    start_R()

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    if not mnemonic_file:
        mnemonic_file = ""  # None cannot be passed to R

    r["save.dinuc"](dirpath, namebase, tabfilename,
                    mnemonic_file=mnemonic_file,
                    clobber=clobber)

def save_html(dirpath, clobber=False, mnemonicfile=None):
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  tables={"":NAMEBASE}, module=MODULE,
                  mnemonicfile=mnemonicfile, title=HTML_TITLE)

## Package entry point
def validate(bedfilename, genomedatadir, dirpath, clobber=False, quick=False,
             replot=False, noplot=False, mnemonic_file=None):
    setup_directory(dirpath)
    if not replot:
        segmentation = Segmentation(bedfilename)
        labels = segmentation.labels

        with Genome(genomedatadir) as genome:
            nuc_counts, dinuc_counts = \
                        calc_nucleotide_frequencies(segmentation,
                                                    genome, quick=quick)

        save_tab(labels, nuc_counts, dinuc_counts, dirpath, clobber=clobber)

    if not noplot:
        save_plot(dirpath, clobber=clobber, mnemonic_file=mnemonic_file)

    save_html(dirpath, clobber=clobber, mnemonicfile=mnemonic_file)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] BEDFILE GENOMEDATADIR"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("--clobber", action="store_true",
                      dest="clobber", default=False,
                      help="Overwrite existing output files if the specified"
                      " directory already exists.")
    parser.add_option("--quick", action="store_true",
                      dest="quick", default=False,
                      help="Compute values only for one chromosome.")
    parser.add_option("--replot", action="store_true",
                      dest="replot", default=False,
                      help="Load data from output tab files and"
                      " regenerate plots instead of recomputing data")
    parser.add_option("--noplot", action="store_true",
                      dest="noplot", default=False,
                      help="Do not generate plots")

    parser.add_option("--mnemonic-file", dest="mnemonic_file",
                      default=None,
                      help="If specified, labels will be shown using"
                      " mnemonics found in this file")
    parser.add_option("-o", "--outdir",
                      dest="outdir", default="%s" % MODULE,
                      help="File output directory (will be created"
                      " if it does not exist) [default: %default]")

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

    validate(bedfilename, genomedatadir, options.outdir,
             clobber=options.clobber, quick=options.quick,
             replot=options.replot, noplot=options.noplot,
             mnemonic_file=options.mnemonic_file)

if __name__ == "__main__":
    sys.exit(main())
