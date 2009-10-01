#!/usr/bin/env python
from __future__ import division, with_statement

"""
validate_all: Evaluate the characteristics of a given
segmentation with respect to a predefined set of genomic landmarks.

This tool, at least in its first incarnation, is intended to provide
spot-checking functionality, giving an overview of the characteristics
of a particular segmentation. The tool performs a series of different
types of analyses, summarizing the results in an HTML report. By
default, all validations are performed, though a subset may be
selected using the command line options listed below.

OUTDIR: A unique name associated with this evaluation. All of the results
are stored in a directory named <root>, with output filenames as
described below.

BEDFILE: The name of a file that contains a segmentation. The file must be
in BED format, and may be gzipped. The first three columns specify
genomic coordinates and the fourth column indicates the segment name.
This name is interpreted as a "state" identifier. Subsequent columns
in the segmentation BED file are ignored.

GENOMEDATADIR: A directory that contains a number of genomedata files

XXX: assumes that segmentation is non-overlapping
"""

__version__ = "$Revision$"

# Copyright 2008 Michael M. Hoffman <mmh1@washington.edu>

import imp
import os
import sys

from contextlib import contextmanager
from inspect import getargspec, isfunction

from .common import die, PKG, setup_directory, tab2html

# If you add a new module to segtools, add the MODULE name
#  (the name of the directory and the field of the MODULE global
#  in that module's validate.py script) to the list below.
MODULES = ["length_distribution",
           "signal_distribution",
           "feature_aggregation",
           "nucleotide_frequency",
           "tss_overlap",
           "label_transition"]

NAMEBASE_SUMMARY = "summary"


########################## HTML FUNCTIONS #########################

def write_html_summary(dirpath, namebase, bedfilename, genomedatadirname, filenames, segmentation):
    with html_writer(dirpath, namebase, bedfilename, genomedatadirname, 
                     filenames, segmentation) as writer:
        write_html_length(dirpath, writer)
        write_html_distributions(dirpath, writer)
        write_html_nucleotide(dirpath, writer)    
        write_html_tss(dirpath, writer)                     
        write_html_aggregation(dirpath, writer)

## A generator for the html writing to allow for 
## selectively writing specific html function output
@contextmanager
def html_writer(dirpath, namebase, bedfilename, genomedatadirname, filenames, segmentation):
    outfilename = make_filename(dirpath, namebase, EXT_HTML)
    with open(outfilename, "w") as outfile:
        write_html_header(outfile, bedfilename, genomedatadirname, filenames, segmentation)
        yield outfile
        write_html_footer(outfile)

    return

########################## MAIN FUNCTIONS #########################

## Main function to handle complete script
## Writes complete HTML summary after running all validators
def validate_all(outdirname, validators={}, arguments={}):
    """
    outdirname: folder to output all results to

    validators: dict
        key: validator module name
        value: validate function to call with arguments

    arguments: dict
        key: argument name (must case-sensitive match function argument names)
        value: argument value
        
    The following arguments (exact case required) will automatically
    be supplied to each validator, if required:
        bedfilename
        genomedatadir
        tssfilename
        dirpath
        clobber
        quick
        calc_ranges
        group_labels
    """
    setup_directory(outdirname)
    
    for func_name, func in validators.iteritems():
	assert isfunction(func)
        print "Running %s" % func_name

        # Set output directory
        arguments["dirpath"] = os.path.join(outdirname, func_name)

        # Assemble func argument dict
        func_args = {}
        args, varargs, varkw, defaults = getargspec(func)
        for arg in args:
            # XXX: If the arg has a default value, don't require a value
            #   to be passed in
            try:
                func_args[arg] = arguments[arg]
                # XXX: Not sure if this should die or silently pass the None
                if arguments[arg] is None:
                    raise KeyError
            except KeyError:
                die(("Module: %s expected argument: %s, but argument"+
                    " not specified") % (func_name, arg))

        # All arguments are satisfied; call function!
        func(**func_args)
        
    #write_html_summary(dirpath, NAMEBASE_SUMMARY, bedfilename, genomedatadirname, filenames, segmentation)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] OUTDIR"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Module options")
    group.add_option("-m", "--module", action="append", 
                     dest="module", default=[],
                     help="Run this segtools validation module on the"
                     " provided data. Multiple modules can be run."
                     " Specify -m or --module before each module name.")
    group.add_option("--bedfile", dest="bedfilename",
                     type="string", default=None,
                     help="BED-formatted segmentation file.")
    group.add_option("--genome-data-dir", dest="genomedatadir",
                     type="string", default=None,
                     help="Directory containing hdf5 chromosome data files.")
    group.add_option("--tss-file", dest="tssfilename",
                     type="string", default=None,
                     help="TSS file in GFF format for the TSS/overlap"
                     " analyses.")
    group.add_option("--clobber", action="store_true",
                     dest="clobber", default=False,
                     help="Overwrite existing output files if the specified"
                     " directory already exists. [default: %default]")
    group.add_option("--quick", action="store_true", 
                     dest="quick", default=False,
                     help="Compute values only for one chromosome."
                     " [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "signal_distribution options")
    group.add_option("--group-labels", action="store_true", 
                      dest="group_labels", default=False,
                      help="Group track distributions over all labels."
                      " BEDFILE will be ignored")
    group.add_option("--calc-ranges", action="store_true", 
                      dest="calc_ranges", default=False,
                      help="Calculate ranges for distribution plots from"
                      " segmentation data (slower) instead of using whole"
                      " genome data (default).")
    parser.add_option_group(group)
    

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Invalid number of arguments")

    if len(options.module) == 0:
        parser.error("Expected at least one module to run. Examples: " +
                     "length_distribution, signal_distribution, " +
                     "nucleotide_frequency, etc.")
        
    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    outdirname = args[0]

    # Save all the options in a normal dict form to make it easier
    # to retrieve the values for the validator
    arguments = {}
    arguments["bedfilename"] = options.bedfilename
    arguments["genomedatadir"] = options.genomedatadir
    arguments["tssfilename"] = options.tssfilename
    arguments["clobber"] = options.clobber
    arguments["quick"] = options.quick
    arguments["group_labels"] = options.group_labels
    arguments["calc_ranges"] = options.calc_ranges
    
    validators = {}
    for modulename in options.module:
        if modulename in MODULES:
            # XXX: Make this less ugly
            fullname = "%s.%s.validate" % (PKG, modulename)
            try:
                __import__(fullname)
                module = sys.modules[fullname]
            except ImportError:
                die("Failed to load module: %s" % fullname)

            # Store module name and validate function
            validators[modulename] = module.validate
        else:
            die("Could not find %s in MODULES list. Make sure to add new" +
                " new modules to the MODULES list in validate_all.py.")

    validate_all(outdirname, validators, arguments)
    
if __name__ == "__main__":
    sys.exit(main())
