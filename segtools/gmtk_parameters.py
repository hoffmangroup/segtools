#!/usr/bin/env python

"""
gmtk_parameters.py
"""

from __future__ import division, with_statement

__version__ = "$Revision$"

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE = "gmtk_parameters"

import os
import sys

from numpy import array

from . import log
from .common import die, r_call, r_source, setup_directory
from .html import save_html_div
from .transition import save_plot, save_graph
from .mnemonics import create_mnemonic_file
from .signal_distribution import save_stats_plot

NAMEBASE = "%s" % MODULE
NAMEBASE_GRAPH = os.extsep.join([NAMEBASE, "graph"])
NAMEBASE_STATS = os.extsep.join([NAMEBASE, "stats"])

HTML_TEMPLATE_FILENAME = "gmtk_div.tmpl"
HTML_TITLE = "GMTK Theoretical Parameters"

P_THRESH = 0.15  # Default
Q_THRESH = 0.0

def start_R():
    r_source("common.R")
    r_source("transition.R")

def get_default_labels(num_labels):
    """Generate default labels (0 through num_labels-1)"""
    return dict([(val, str(val)) for val in range(0, num_labels)])

def load_gmtk_transitions(gmtk_file):
    """Loads probabilites from a gmtk_file, through R.

    Returns probs as a numpy.array
    """
    start_R()

    r_data = r_call("read.gmtk.transition", gmtk_file)
    # Rpy automatically transposes, so need to transpose it back
    probs = array(r_data, dtype="double").transpose()
    num_labels = probs.shape[0]
    labels = get_default_labels(num_labels)

    return labels, probs

def save_html(dirpath, gmtk_file, p_thresh, q_thresh, clobber=False):
    extra_namebases = {"graph": NAMEBASE_GRAPH,
                       "stats": NAMEBASE_STATS}

    if p_thresh > 0:
        thresh = "P > %s" % p_thresh
    elif q_thresh > 0:
        thresh = "P above %.2th quantile" % q_thresh

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, extra_namebases=extra_namebases,
                  title=HTML_TITLE, thresh=thresh,
                  gmtk_file=os.path.basename(gmtk_file))

## Package entry point
def validate(gmtk_file, dirpath, p_thresh=P_THRESH, q_thresh=Q_THRESH,
             noplot=False, nograph=False, gene_graph=False, clobber=False,
             translation_file=None, allow_regex=False, mnemonic_file=None,
             create_mnemonics=False, verbose=True):

    setup_directory(dirpath)

    if not os.path.isfile(gmtk_file):
        die("Could not find GMTK file: %s" % gmtk_file)

    # XXX: why does this print "getbuffer" a few times? debug R
    # read.gmtk.transition() code
    log("Loading gmtk transitions...", verbose, end="")

    labels, probs = load_gmtk_transitions(gmtk_file)
    log(" done", verbose)

    # If mnemonics weren't specified, let's create a mnemonic file
    if mnemonic_file is None and create_mnemonics:
        mnemonic_file = \
            create_mnemonic_file(gmtk_file, dirpath, verbose=verbose,
                                 clobber=clobber, gmtk=True)

    if not noplot:
        save_plot(dirpath, namebase=NAMEBASE, filename=gmtk_file,
                  verbose=verbose, clobber=clobber, gmtk=True,
                  mnemonic_file=mnemonic_file)
        save_stats_plot(dirpath, namebase=NAMEBASE_STATS, filename=gmtk_file,
                        clobber=clobber, gmtk=True,
                        mnemonic_file=mnemonic_file,
                        translation_file=translation_file,
                        allow_regex=allow_regex, verbose=verbose)

    if not nograph:
        save_graph(labels, probs, dirpath, clobber=clobber,
                   p_thresh=p_thresh, q_thresh=q_thresh,
                   gene_graph=gene_graph, mnemonic_file=mnemonic_file,
                   namebase=NAMEBASE_GRAPH, verbose=verbose)

    save_html(dirpath, gmtk_file, p_thresh=p_thresh, q_thresh=q_thresh,
              clobber=clobber)


def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] PARAMSFILE"
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
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate transition plots")
    group.add_option("--nograph", action="store_true",
                     dest="nograph", default=False,
                     help="Do not generate transition graph")
    group.add_option("--create-mnemonics", action="store_true",
                     dest="create_mnemonics", default=False,
                     help="If mnemonics are not specified, they will be"
                     " created and used for plotting")
    group.add_option("--allow-regex", dest="allow_regex",
                     default=False, action="store_true",
                     help="If a file is specified with --trackname-translation"
                     ", the mappings in this file will be interpreted as"
                     " regular expressions instead of exact mappings."
                     " Thus, all underscores could be converted to periods"
                     " with the single line: `_<TAB>.`. Mappings will be"
                     " applied sequentially, so a second line of:"
                     " `\\.<TAB>#` would (together with the first line)"
                     " convert all periods and underscores to pounds.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    group.add_option("-m", "--mnemonic-file", dest="mnemonic_file",
                     default=None, metavar="FILE",
                     help="If specified, labels will be shown using"
                     " mnemonics found in this file")
    group.add_option("-t", "--trackname-translation", dest="translation_file",
                     default=None, metavar="FILE",
                     help="Should be a file with rows <old-trackname>"
                     "<TAB><new-trackname>. Tracknames will be translated"
                     " using this mapping before plotting the stats plot."
                     " By default, <old-trackname> must exactly match"
                     " the name of a track, but --allow-regex provides"
                     " more flexibility.")
    group.add_option("-o", "--outdir",
                     dest="outdir", default="%s" % MODULE, metavar="DIR",
                     help="File output directory (will be created"
                     " if it does not exist) [default: %default]")
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
    group.add_option("--gene-graph", dest="gene_graph",
                     default=False, action="store_true",
                     help="Make each node of the graph a reference to a .ps"
                     " image an \"image\" subdirectory. Currently, these .ps"
                     " files need to be made separately.")
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

    args = [args[0], options.outdir]
    kwargs = {"p_thresh": options.p_thresh,
              "q_thresh": options.q_thresh,
              "clobber": options.clobber,
              "verbose": options.verbose,
              "noplot": options.noplot,
              "nograph": options.nograph,
              "create_mnemonics": options.create_mnemonics,
              "gene_graph": options.gene_graph,
              "translation_file": options.translation_file,
              "allow_regex": options.allow_regex,
              "mnemonic_file": options.mnemonic_file}
    validate(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
