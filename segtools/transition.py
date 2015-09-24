#!/usr/bin/env python
from __future__ import division, with_statement

"""
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
MODULE = "transition"


import os
import sys


from numpy import empty, loadtxt, where, zeros
from subprocess import call

from . import log, Segmentation, die, RInterface, open_transcript, \
     add_common_options
from .common import check_clobber, get_ordered_labels, make_dotfilename, \
     make_pdffilename, make_pngfilename, make_namebase_summary, \
     make_tabfilename, map_mnemonics, setup_directory, tab_saver
from .html import save_html_div
from .mnemonics import create_mnemonic_file
from .version import __version__

NAMEBASE = "%s" % MODULE
NAMEBASE_SUMMARY = make_namebase_summary(NAMEBASE)
NAMEBASE_GRAPH = os.extsep.join([NAMEBASE, "graph"])

HTML_TEMPLATE_FILENAME = "transition_div.tmpl"
HTML_TITLE = "Segment label transitions"

P_THRESH = 0.15  # Default
Q_THRESH = 0.0

BASE_WEIGHT = 0.3  # Base edge/arrow weight
EDGE_WEIGHT = 2.7
ARROW_WEIGHT = 1.7

R = RInterface(["common.R", "transition.R"])

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
                if row1['end'] == row2['start']:
                    label_key1 = row1['key']
                    label_key2 = row2['key']
                    counts[label_key1, label_key2] += 1
                    if label_key1 == label_key2:
                        if self_transition_count == 0:
                            log("WARNING: unexpected self-transition"
                                "%s -> %s\n\tRow 1: %s\n\tRow 2: %s" % \
                                    (labels[label_key1], labels[label_key2],
                                     row1, row2))
                        self_transition_count += 1

            row1 = row2  # Inch forward

    if self_transition_count > 0:
        log("WARNING: %d self-transitions observed" % self_transition_count)

    # Row-normalize to turn counts into probabilities
    probs = empty(counts.shape, dtype='double')
    for i, row in enumerate(counts):
        probs[i] = row / sum(row)

    return labels, probs

def load_transition_probs(dirpath, namebase=NAMEBASE):
    filename = make_tabfilename(dirpath, namebase)
    return loadtxt(filename)

def save_tab(labels, probs, dirpath, clobber=False, verbose=True):
    ordered_keys, labels = get_ordered_labels(labels)

    # Get fieldnames in order
    fieldnames = list(labels[key] for key in ordered_keys)
    with tab_saver(dirpath, NAMEBASE, fieldnames, verbose=verbose,
                   clobber=clobber, header=False) as saver:
        for start_key in ordered_keys:
            prob_row = probs[start_key]
            row = {}
            for end_key in ordered_keys:
                row[labels[end_key]] = "%.5f" % prob_row[end_key]
            saver.writerow(row)

def save_html(dirpath, p_thresh, q_thresh, clobber=False, verbose=True):
    extra_namebases = {"graph": NAMEBASE_GRAPH}

    if p_thresh > 0:
        thresh = "P >= %s" % p_thresh
    elif q_thresh > 0:
        thresh = "P above %.2th quantile" % q_thresh

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, extra_namebases=extra_namebases,
                  title=HTML_TITLE, thresh=thresh, verbose=verbose)

def save_plot(dirpath, namebase=NAMEBASE, filename=None, ddgram=False,
              clobber=False, mnemonic_file=None, verbose=False, gmtk=False,
              transcriptfile=None, ropts=None):
    """
    if filename is specified, it is used instead of dirpath/namebase.tab
    """
    R.start(verbose=verbose, transcriptfile=transcriptfile)

    if filename is None:
        filename = make_tabfilename(dirpath, namebase)

    if not os.path.isfile(filename):
        die("Unable to find tab file: %s" % filename)

    R.plot("save.transition", dirpath, namebase, filename, clobber=clobber,
           mnemonic_file=mnemonic_file, ddgram=ddgram,
           gmtk=gmtk, ropts=ropts)

def save_graph(labels, probs, dirpath, q_thresh=Q_THRESH, p_thresh=P_THRESH,
               clobber=False, mnemonic_file=None, fontname="Helvetica",
               lenient_thresh=True, gene_graph=False, namebase=NAMEBASE_GRAPH,
               verbose=True):
    assert labels is not None and probs is not None
    try:
        import pygraphviz as pgv
    except ImportError:
        log("Unable to load PyGraphviz library. Skipping transition graph")
        return

    dotfilename = make_dotfilename(dirpath, namebase)
    check_clobber(dotfilename, clobber)
    pngfilename = make_pngfilename(dirpath, namebase)
    check_clobber(pngfilename, clobber)
    pdffilename = make_pdffilename(dirpath, namebase)
    check_clobber(pdffilename, clobber)

    # Replace labels with mnemonic labels, if mnemonics are given
    mnemonics = map_mnemonics(labels, mnemonic_file)
    ordered_keys, labels = get_ordered_labels(labels, mnemonics)

    # Threshold
    if q_thresh > 0:
        quantile = R.call('matrix.find_quantile', probs, q_thresh)[0]
        log("Removing edges below %.4f" % float(quantile), verbose)

        probs[probs < quantile] = 0
    elif p_thresh > 0:
        log("Removing connections below %.4f" % p_thresh, verbose)

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

    G.node_attr.update(fontname=fontname)
    if gene_graph:
        G.node_attr.update(shape="plaintext")
        ps_dir = os.path.join(os.path.split(dirpath)[0], "images")
        for node in G.nodes():
            node.attr["image"] = "%s/%s.ps" % (ps_dir, str(node))
            node.attr["label"] = " "

    G.write(dotfilename)

    log("Drawing graphs...", verbose, end="")

    G.layout()

    try:
        G.draw(pngfilename)
    except:
        log("Error: Failed to draw png graph")

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
        log("Error: Failed to draw pdf graph: %s" % str(e))

    log(" done", verbose)

## Package entry point
def validate(bedfilename, dirpath, ddgram=False, p_thresh=P_THRESH,
             q_thresh=Q_THRESH, replot=False, noplot=False, nograph=False,
             gene_graph=False, clobber=False, gmtk=False, verbose=True,
             mnemonic_file=None, ropts=None):
    setup_directory(dirpath)

    if gmtk:
        from .gmtk_parameters import load_gmtk_transitions
        labels, probs = load_gmtk_transitions(bedfilename)
        if mnemonic_file is None:
            mnemonic_file = create_mnemonic_file(bedfilename, dirpath,
                                                 clobber=clobber,
                                                 gmtk=gmtk, verbose=verbose)
    else:
        segmentation = Segmentation(bedfilename, verbose=verbose)
        assert segmentation is not None

        if replot:
            probs = load_transition_probs(dirpath)
            labels = segmentation.labels
        else:
            # Calculate transition probabilities for each label
            labels, probs = calc_transitions(segmentation)


    if not replot:
        save_tab(labels, probs, dirpath, clobber=clobber, verbose=verbose)

    if not noplot:
        with open_transcript(dirpath, MODULE) as transcriptfile:
            save_plot(dirpath, ddgram=ddgram, verbose=verbose,
                      clobber=clobber, mnemonic_file=mnemonic_file,
                      transcriptfile=transcriptfile, ropts=ropts)

    if not nograph:
        save_graph(labels, probs, dirpath, clobber=clobber,
                   p_thresh=p_thresh, q_thresh=q_thresh,
                   gene_graph=gene_graph, mnemonic_file=mnemonic_file,
                   verbose=verbose)

    save_html(dirpath, p_thresh=p_thresh, q_thresh=q_thresh,
              clobber=clobber, verbose=verbose)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] SEGMENTATION"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'replot', 'noplot'])
    group.add_option("--nograph", action="store_true",
                     dest="nograph", default=False,
                     help="Do not generate transition graph")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Transition frequency plot options")
    group.add_option("--dd", "--dendrogram", dest="ddgram",
                     default=False, action="store_true",
                     help="include dendrogram along edge of levelplot"
                     " [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Transition graph options")
    group.add_option("-P", "--prob-threshold", dest="p_thresh",
                     type="float", default=P_THRESH, metavar="VAL",
                     help="ignore all transitions with probabilities below"
                     " this absolute threshold [default: %default]")
    group.add_option("-Q", "--quantile-threshold", dest="q_thresh",
                     type="float", default=Q_THRESH, metavar="VAL",
                     help="ignore transitions with probabilities below this"
                     " probability quantile [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Non-segmentation files")
    group.add_option("-g", "--gmtk", action="store_true",
                     dest="gmtk", default=False,
                     help="The SEGMENTATION argument will instead be treated"
                     " as a GMTK parameter file. If a mnemonic file is not"
                     " specified, one will be created and used.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

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

    kwargs = dict(options.__dict__)
    args = [args[0], kwargs.pop('outdir')]
    validate(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
