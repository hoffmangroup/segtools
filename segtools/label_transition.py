#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
label_transition.py:

Provides command-line and package entry points for analyzing the observed
segmentation label transitions in the given BED-formatted segmentation.


Accepts an input file containing a matrix of transition
probabilities and generates a heatmap of the matrix
in pdf and png formats.

Optionally generates a heatmap of the ratio of reciprocal
transition probabilities (a->b/b->a)

Optionally generates a graphviz image of a thresholded
form of the transition matrix

XXX: Graphiz currently assumes integer labels in sorted order in file

"""

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE = "label_transition"


import os
import sys


from contextlib import contextmanager
from functools import partial
from numpy import array, empty, sort, where, zeros
from subprocess import call
from rpy2.robjects import r, rinterface, numpy2ri

# XXX: Do this without the kludgy constants
from .common import check_clobber, die, get_ordered_labels, image_saver, load_segmentation, load_genome, make_dotfilename, make_filename, make_pdffilename, make_pngfilename, make_namebase_summary, make_tabfilename, map_mnemonics, OutputMasker, r_source, SEGMENT_START_COL, SEGMENT_END_COL, SEGMENT_LABEL_KEY_COL, setup_directory, tab_saver

from .html import save_html_div


NAMEBASE = "%s" % MODULE
NAMEBASE_SUMMARY = make_namebase_summary(NAMEBASE)
NAMEBASE_GRAPH = os.extsep.join([NAMEBASE, "graph"])
NAMEBASE_PARAMS = os.extsep.join([NAMEBASE, "params"])

HTML_TEMPLATE_FILENAME = "transition_div.tmpl"
HTML_TITLE = "Segment label transitions"

P_THRESH = 0.15  # Default
Q_THRESH = 0.0

BASE_WEIGHT = 0.3  # Base edge/arrow weight
EDGE_WEIGHT = 2.7
ARROW_WEIGHT = 1.7

PNG_SIZE = 600  # px
PNG_SIZE_PARAMS = 800  # px

def start_R():
    r_source("common.R")
    r_source("transition.R")

## Returns array of row-normalized probabilities corresponding
## to transition probabilities between labels
## Start labels are rows, end labels are cols
def calc_transitions(segmentation):
    labels = segmentation.labels
    self_transition_count = 0
    # NxN matrix for transition counts
    counts = zeros((len(labels), len(labels)), dtype="int")
    for segments in segmentation.chromosomes.itervalues():
        # Inch along rows, looking at all pairs
        row_it = iter(segments)
        row1 = None
        for row2 in row_it:
            if row1 is not None:  # Skip first loop
                # Increment transition only if directly adjacent
                if row1[SEGMENT_END_COL] == row2[SEGMENT_START_COL]:
                    label_key1 = row1[SEGMENT_LABEL_KEY_COL]
                    label_key2 = row2[SEGMENT_LABEL_KEY_COL]
                    counts[label_key1, label_key2] += 1
                    if label_key1 == label_key2:
                        if self_transition_count == 0:
                            print >>sys.stderr, \
                                ("WARNING: unexpected self-transition"
                                 "%s -> %s\n\tRow 1: %s\n\tRow 2: %s" % \
                                (labels[label_key1], labels[label_key2],
                                 row1, row2))
                        self_transition_count += 1

            row1 = row2  # Inch forward

    if self_transition_count > 0:
        print >>sys.stderr, "WARNING: %d self-transitions observed" % \
            self_transition_count

    # Row-normalize to turn counts into probabilities
    probs = empty(counts.shape, dtype='double')
    for i, row in enumerate(counts):
        probs[i] = row / sum(row)

    return labels, probs


def save_tab(labels, probs, dirpath, clobber=False):
    ordered_keys, labels = get_ordered_labels(labels)

    # Get fieldnames in order
    fieldnames = list(labels[key] for key in ordered_keys)
    with tab_saver(dirpath, NAMEBASE, fieldnames,
                   clobber=clobber, header=False) as saver:
        for start_key in ordered_keys:
            prob_row = probs[start_key]
            row = {}
            for end_key in ordered_keys:
                row[labels[end_key]] = "%.5f" % prob_row[end_key]
            saver.writerow(row)

def save_summary_tab(labels, probs, dirpath, mnemonics=None, clobber=False):
    ordered_keys, labels = get_ordered_labels(labels, mnemonics)
    # Get fieldnames in order
    assert "label" not in labels.keys()
    fieldnames = ["label"] + list(labels[key] for key in ordered_keys)
    with tab_saver(dirpath, NAMEBASE_SUMMARY, fieldnames,
                   clobber=clobber) as saver:
        for start_key in ordered_keys:
            prob_row = probs[start_key]
            row = {}
            row["label"] = labels[start_key]
            for end_key in ordered_keys:
                row[labels[end_key]] = "%.5f" % prob_row[end_key]
            saver.writerow(row)

def save_html(dirpath, p_thresh, q_thresh, clobber=False):
    extra_namebases = {"graph": NAMEBASE_GRAPH}

    if p_thresh > 0:
        thresh = "P > %s" % p_thresh
    elif q_thresh > 0:
        thresh = "P above %.2th quantile" % q_thresh

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, extra_namebases=extra_namebases,
                  title=HTML_TITLE, thresh=thresh)

def save_plot(dirpath, ddgram=False, clobber=False, mnemonics=[],
              verbose=False):
    start_R()

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    with image_saver(dirpath, NAMEBASE,
                     height=PNG_SIZE,
                     width=PNG_SIZE,
                     clobber=clobber,
                     verbose=verbose):
        r.plot(r["plot.transition"](tabfilename,
                                    mnemonics=mnemonics,
                                    ddgram=ddgram))

def plot_gmtk_transitions(filename, dirpath, clobber=False, verbose=False,
                          mnemonics=[]):
    """Plots transition probabilities from gmtk params file.
    """
    assert filename is not None and os.path.isfile(filename)

    start_R()

    # Plot transition probabilies
    with image_saver(dirpath, NAMEBASE,
                     height=PNG_SIZE,
                     width=PNG_SIZE,
                     clobber=clobber,
                     verbose=verbose):
        r.plot(r["plot.gmtk.transition"](filename,
                                         mnemonics=mnemonics))

    # Plot track parameters
    with image_saver(dirpath, NAMEBASE_PARAMS,
                     height=PNG_SIZE_PARAMS,
                     width=PNG_SIZE_PARAMS,
                     clobber=clobber,
                     verbose=verbose):
        r.plot(r["plot.gmtk.params"](filename,
                                     mnemonics=mnemonics))

def save_graph(labels, probs, dirpath, q_thresh=Q_THRESH, p_thresh=P_THRESH,
               clobber=False, mnemonics=[], fontname="Helvetica",
               lenient_thresh=True, gene_graph=False):
    assert labels is not None and probs is not None
    try:
        import pygraphviz as pgv
    except ImportError:
        print >>sys.stderr, ("Unable to load PyGraphviz library."
                             " Skipping transition graph")
        return

    dotfilename = make_dotfilename(dirpath, NAMEBASE_GRAPH)
    check_clobber(dotfilename, clobber)
    pngfilename = make_pngfilename(dirpath, NAMEBASE_GRAPH)
    check_clobber(pngfilename, clobber)
    pdffilename = make_pdffilename(dirpath, NAMEBASE_GRAPH)
    check_clobber(pdffilename, clobber)

    # Replace labels with mnemonic labels, if mnemonics are given
    ordered_keys, labels = get_ordered_labels(labels, mnemonics)

    # Threshold
    if q_thresh > 0:
        quantile = r['matrix.find_quantile'](probs, q_thresh)[0]
        print >>sys.stderr, "Removing edges below %.4f" % float(quantile)
        probs[probs < quantile] = 0
    elif p_thresh > 0:
        print >>sys.stderr, "Removing connections below %.4f" % p_thresh
        if lenient_thresh:
            row_wise_remove = (probs < p_thresh)
            col_wise = probs / probs.sum(axis=0) # Col-normalize
            col_wise_remove = (col_wise < p_thresh)
            probs[row_wise_remove & col_wise_remove] = 0
        else:
            probs[probs < p_thresh] = 0

    # Create graph out of non-zero edges
    G = pgv.AGraph(strict=False, directed=True)

    rows, cols = where(probs > 0)

    max_val = probs.max()
    min_val = probs[probs > 0].min()
    weights = (probs - min_val)/(max_val-min_val) + BASE_WEIGHT
    weights[weights < BASE_WEIGHT] = 0
    for row, col in zip(rows, cols):
        weight = weights[row, col]
        G.add_edge(labels[row], labels[col],
                   penwidth=str(EDGE_WEIGHT * weight),
                   arrowsize=str(ARROW_WEIGHT * weight))

    if gene_graph:
        G.node_attr.update(fontname=fontname, shape="plaintext")
        ps_dir = os.path.join(os.path.split(dirpath)[0], "images")
        for node in G.nodes():
            node.attr["image"] = "%s/%s.ps" % (ps_dir, str(node))
            node.attr["label"] = " "

    G.write(dotfilename)

    print >>sys.stderr, "Drawing graphs...",
    G.layout()

    try:
        G.draw(pngfilename)
    except:
        print >>sys.stderr, "Failed to draw png graph"

    try:
        # Try to generate pdf from dot
        if gene_graph:
            layout_prog = "dot"
        else:
            layout_prog = "neato"

        cmd = " ".join([layout_prog, "-Tps2", dotfilename, "|",
                        "ps2pdf", "-dAutoRotatePages=/None", "-", pdffilename])
        code = call(cmd, shell=True)
        if code != 0:
            raise Exception()
    except Exception, e:
        print >>sys.stderr, "Failed to draw pdf graph: %s" % str(e)

    print >>sys.stderr, "done"

def create_mnemonic_file(gmtk_file, dirpath, clobber=False, namebase=None):
    """Generate a mnemonic file for GMTK params with R clustering code.

    Calls R code that writes mnemonics to a file.
    Returns name of created mnemonic file
    """
    start_R()

    if namebase is None:
        gmtk_filebase = os.path.basename(gmtk_file)
        if gmtk_filebase.endswith(".params"):
            namebase = gmtk_filebase[:-7]
        else:
            namebase = gmtk_filebase

    mnemonic_base = os.extsep.join([namebase, "mnemonics"])

    filename = make_tabfilename(dirpath, mnemonic_base)
    check_clobber(filename, clobber)

    # Create mnemonic file via R
    try:
        r_filename = r["make.gmtk.mnemonic.file"](gmtk_file, filename=filename)
        mnemonicfilename = str(r_filename[0])

        print >>sys.stderr, "Generated mnemonic file: %s" % mnemonicfilename
        return mnemonicfilename
    except rinterface.RRuntimeError:
        print >>sys.stderr, ("ERROR: Failed to create mnemonic file."
                             " Continuing without mnemonics.")
        return None

def get_num_gmtk_labels(gmtk_file):
    """Load number of labels from gmtk params file."""
    start_R()

    num_labels = int(r["num.gmtk.labels"](gmtk_file)[0])
    return num_labels

def get_default_labels(num_labels):
    """Generate default labels (0 through num_labels-1)"""
    return dict([(val, str(val)) for val in range(0, num_labels)])

def load_gmtk_transitions(gmtk_file):
    """Loads probabilites from a gmtk_file, through R.

    Returns probs as a numpy.array
    """
    start_R()

    r_data = r["read.gmtk.transition"](gmtk_file)
    # Rpy automatically transposes, so need to transpose it back
    probs = array(r_data, dtype="double").transpose()
    return probs


## Package entry point
def validate_gmtk(gmtk_file, dirpath, ddgram=False, p_thresh=P_THRESH,
                  q_thresh=Q_THRESH, noplot=False, nograph=False,
                  gene_graph=False, clobber=False, mnemonicfilename=None):
    setup_directory(dirpath)

    if not os.path.isfile(gmtk_file):
        die("Could not find GMTK file: %s" % gmtk_file)

    # If mnemonics weren't specified, let's create a mnemonic file!
    if mnemonicfilename is None:
        mnemonicfilename = create_mnemonic_file(gmtk_file, dirpath,
                                                clobber=clobber)

    num_labels = get_num_gmtk_labels(gmtk_file)
    labels = get_default_labels(num_labels)
    mnemonics = map_mnemonics(labels, mnemonicfilename)

    if not noplot:
        plot_gmtk_transitions(gmtk_file, dirpath, mnemonics=mnemonics,
                              clobber=clobber)

    if not nograph:
        probs = load_gmtk_transitions(gmtk_file)
        save_graph(labels, probs, dirpath, clobber=clobber,
                   p_thresh=p_thresh, q_thresh=q_thresh,
                   gene_graph=gene_graph, mnemonics=mnemonics)

    save_html(dirpath, p_thresh=p_thresh, q_thresh=q_thresh,
              clobber=clobber)

## Package entry point
def validate(bedfilename, dirpath, ddgram=False, p_thresh=P_THRESH,
             q_thresh=Q_THRESH, noplot=False, nograph=False,
             gene_graph=False, clobber=False, mnemonicfilename=None):
    setup_directory(dirpath)

    segmentation = load_segmentation(bedfilename)
    assert segmentation is not None

    # Calculate transition probabilities for each label
    labels, probs = calc_transitions(segmentation)

    mnemonics = map_mnemonics(labels, mnemonicfilename)
    save_summary_tab(labels, probs, dirpath,
                     clobber=clobber, mnemonics=mnemonics)
    save_tab(labels, probs, dirpath, clobber=clobber)
    if not noplot:
        save_plot(dirpath, ddgram=ddgram,
                  clobber=clobber, mnemonics=mnemonics)

    if not nograph:
        save_graph(labels, probs, dirpath, clobber=clobber,
                   p_thresh=p_thresh, q_thresh=q_thresh,
                   gene_graph=gene_graph, mnemonics=mnemonics)

    save_html(dirpath, p_thresh=p_thresh, q_thresh=q_thresh,
              clobber=clobber)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] BEDFILE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("--clobber", action="store_true",
                      dest="clobber", default=False,
                      help="Overwrite existing output files if the specified"
                      " directory already exists.")
    parser.add_option("--noplot", action="store_true",
                      dest="noplot", default=False,
                      help="Do not generate transition plots")
    parser.add_option("--nograph", action="store_true",
                      dest="nograph", default=False,
                      help="Do not generate transition graph")
    parser.add_option("--mnemonic-file", dest="mnemonic_file",
                      default=None,
                      help="If specified, labels will be shown using"
                      " mnemonics found in this file")
    parser.add_option("-o", "--outdir",
                      dest="outdir", default="%s" % MODULE,
                      help="File output directory (will be created"
                      " if it does not exist) [default: %default]")

    group = OptionGroup(parser, "Transition frequency plot options")
    group.add_option("--dd", "--dendrogram", dest="ddgram",
                     default=False, action="store_true",
                     help="include dendrogram along edge of levelplot"
                     " [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Transition graph options")
    group.add_option("-p", "--prob-threshold", dest="p_thresh",
                     type="float", default=P_THRESH,
                     help="ignore all transitions with probabilities below"
                     " this absolute threshold [default: %default]")
    group.add_option("-q", "--quantile-threshold", dest="q_thresh",
                     type="float", default=Q_THRESH,
                     help="ignore transitions with probabilities below this"
                     " probability quantile [default: %default]")
    group.add_option("-g", "--gene-graph", dest="gene_graph",
                     default=False, action="store_true",
                     help="Make each node of the graph a reference to a .ps"
                     " image an \"image\" subdirectory. Currently, these .ps"
                     " files need to be made separately.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Non-segmentation files")
    parser.add_option("--gmtk-params", dest="gmtk_file",
                      default=None,
                      help="If specified, parameters in the given GMTK"
                      " file will be used to generate plots instead of the"
                      " observed transitions in the BEDFILE. The BEDFILE"
                      " will not be used")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 1 and options.gmtk_file is None:
        parser.error("Insufficient number of arguments")

    if options.p_thresh > 0 and options.q_thresh > 0:
        parser.error("Cannot specify both absolute and quantile thresholds")
    if options.q_thresh < 0 or options.q_thresh > 1:
        parser.error("Quantile threshold should be in range [0, 1]")
    if options.p_thresh < 0 or options.p_thresh > 1:
        parser.error("Probability threshold should be in range [0, 1]")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)

    if options.gmtk_file is not None:
        validate_gmtk(options.gmtk_file, options.outdir,
                      p_thresh=options.p_thresh, q_thresh=options.q_thresh,
                      clobber=options.clobber, noplot=options.noplot,
                      nograph=options.nograph, gene_graph=options.gene_graph,
                      mnemonicfilename=options.mnemonic_file)

    else:
        bedfilename=args[0]
        validate(bedfilename, options.outdir, ddgram=options.ddgram,
                 p_thresh=options.p_thresh, q_thresh=options.q_thresh,
                 clobber=options.clobber, noplot=options.noplot,
                 nograph=options.nograph, gene_graph=options.gene_graph,
                 mnemonicfilename=options.mnemonic_file)

if __name__ == "__main__":
    sys.exit(main())
