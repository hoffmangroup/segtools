#!/usr/bin/env python
from __future__ import division, with_statement

"""
Provides command-line and package entry points for analyzing nucleotide
and dinucleotide frequencies for each label in the segmentation or annotation
file.

"""

import numpy
import os
import sys
import warnings

from collections import defaultdict
from genomedata import Genome
from numpy import zeros, bincount, array

from . import log, Annotation, die, RInterface, open_transcript, \
     add_common_options, ProgressBar
from .common import make_tabfilename, setup_directory, tab_saver

from .html import save_html_div
from .version import __version__

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

R = RInterface(["common.R", "dinucleotide.R"])

def reduce_nucs(arr):
    # Return array with nucs compressed to 0-4
    return numpy.digitize(arr, [ord(i) for i in 'CGNT'])

def zip_nucs(nuc1, nuc2):
    # Given 2 reduced nuc arrays, return the dinuc value where dinucs are
    # formed from zip(nuc1, nuc2)
    # each reduced nuc array has a range of 0-4
    return nuc1 * 5 + nuc2

## Caclulates nucleotide and dinucleotide frequencies over the specified
## annotation data
def calc_nucleotide_frequencies(annotation, genome,
                                nuc_categories=NUC_CATEGORIES,
                                dinuc_categories=DINUC_CATEGORIES,
                                quick=False, verbose=True):

    # Store categories efficiently as dict
    # from each entry directly to category index
    quick_nuc_categories = defaultdict(dict)
    for index, category in enumerate(nuc_categories):
        for entry in category:
            key = reduce_nucs(array([ord(entry)], dtype='i8'))[0]
            #print("Mapping %s -> %d -> %d" % (entry, key, index))
            assert key not in quick_nuc_categories, \
                   "Error: nucleotide hash for %s is not unique" % entry
            quick_nuc_categories[key] = index
    quick_nuc_categories = dict(quick_nuc_categories)

    quick_dinuc_categories = defaultdict(dict)
    for index, category in enumerate(dinuc_categories):
        for entry in category:
            # Second letter between 0-4, first letter (0-4)*5
            nucs = reduce_nucs(array([ord(i) for i in entry], dtype='i8'))
            key = zip_nucs(nucs[0], nucs[1])
            #print("Mapping %s%s -> %d -> %d" % (entry[0], entry[1], key, index))
            assert key not in quick_dinuc_categories, \
                   "Error: dinucleotide hash for %s%s is not unique" % \
                   (entry[0], entry[1])
            quick_dinuc_categories[key] = index
    quick_dinuc_categories = dict(quick_dinuc_categories)

    # Store counts of each (di)nucleotide observed
    # separated by label
    nuc_counts = defaultdict(dict)
    dinuc_counts = defaultdict(dict)

    labels = annotation.labels
    for label_key in labels.iterkeys():
        nuc_counts[label_key] = zeros(len(nuc_categories)+1, dtype=numpy.long)
        dinuc_counts[label_key] = zeros(len(dinuc_categories)+1,
                                        dtype=numpy.long)

    # Count (di)nucleotides over annotation
    for chromosome in genome:
        chrom = chromosome.name

        try:
            chrom_rows = annotation.chromosomes[chrom]
            # Store entire chromosome's sequence as string in memory for speed
            warnings.simplefilter("ignore")
            log("  %s" % chrom, verbose)
            log("    preparing chromosome sequence...", verbose)
            chrom_sequence = chromosome.seq[0:chromosome.end]
            # Convert to capital letter within integer representation
            chrom_sequence[chrom_sequence >= ord('a')] -= ord('a') - ord('A')
            # Convert to reduced integer representation (0-4)
            chrom_sequence = reduce_nucs(chrom_sequence)
            warnings.resetwarnings()
        except KeyError:
            continue

        log("    counting nucs and dinucs...", verbose)
        if verbose:
            progress = ProgressBar(len(chrom_rows), label="    ")

        for row in chrom_rows:
            if verbose:
                progress.next()

            start = row['start']
            end = row['end']
            label = row['key']

            lab_nuc_counts = nuc_counts[label]
            lab_dinuc_counts = dinuc_counts[label]
            sequence = chrom_sequence[start:end]
            if len(sequence) == 0: continue

            # Count nucs
            val_counts = bincount(sequence)
            for val in val_counts.nonzero()[0]:
                n = val_counts[val]
                # Put in last bin if not found
                index = quick_nuc_categories.get(val, len(lab_nuc_counts) - 1)
                lab_nuc_counts[index] += n

            # Count dinucs efficiently
            if len(sequence) > 1:
                val_counts = bincount(zip_nucs(sequence[:-1], sequence[1:]))
                for val in val_counts.nonzero()[0]:
                    n = val_counts[val]
                    index = quick_dinuc_categories.get(val, len(lab_dinuc_counts) - 1)
                    dinuc_counts[label][index] += n


        if verbose:
            progress.end()
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
             fieldnames=FIELDNAMES, clobber=False, verbose=True):
    with tab_saver(dirpath, NAMEBASE, fieldnames,
                   clobber=clobber, verbose=verbose) as saver:
        for label_key in sorted(labels.keys()):
            row = make_row(labels[label_key],
                           nuc_counts[label_key],
                           dinuc_counts[label_key])
            saver.writerow(row)

def save_plot(dirpath, clobber=False, verbose=True, ropts=None,
              mnemonic_file="", namebase=NAMEBASE, transcriptfile=None):
    R.start(verbose=verbose, transcriptfile=transcriptfile)

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    R.plot("save.dinuc", dirpath, namebase, tabfilename,
           mnemonic_file=mnemonic_file, clobber=clobber, ropts=ropts)

def save_html(dirpath, clobber=False, verbose=True, mnemonicfile=None):
    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  tables={"":NAMEBASE}, module=MODULE, verbose=verbose,
                  mnemonicfile=mnemonicfile, title=HTML_TITLE)

## Package entry point
def validate(filename, genomedatadir, dirpath, clobber=False, quick=False,
             replot=False, noplot=False, mnemonic_file=None, verbose=True,
             ropts=None):
    setup_directory(dirpath)
    if not replot:
        annotation = Annotation(filename, verbose=verbose)
        labels = annotation.labels

        with Genome(genomedatadir) as genome:
            nuc_counts, dinuc_counts = \
                calc_nucleotide_frequencies(annotation, genome,
                                            quick=quick, verbose=verbose)

        save_tab(labels, nuc_counts, dinuc_counts, dirpath,
                 clobber=clobber, verbose=verbose)

    if not noplot:
        with open_transcript(dirpath, MODULE) as transcriptfile:
            save_plot(dirpath, clobber=clobber, verbose=verbose,
                      mnemonic_file=mnemonic_file, ropts=ropts,
                      transcriptfile=transcriptfile)

    save_html(dirpath, clobber=clobber, mnemonicfile=mnemonic_file,
              verbose=verbose)

def parse_options(args):
    from optparse import OptionGroup, OptionParser

    usage = "%prog [OPTIONS] ANNOTATION GENOMEDATAFILE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'quick', 'replot',
                               'noplot'])
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    filename = args[0]
    genomedatadir = args[1]

    kwargs = dict(options.__dict__)
    outdir = kwargs.pop("outdir")
    validate(filename, genomedatadir, outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
