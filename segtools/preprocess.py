#!/usr/bin/env python
"""
If Segtools is taking too long to parse your large segmentation or
annotation file (specify --annotation), this
tool allows you to pre-process that file only once
and have it load much faster in the future.
INFILE will be parsed to create a special binary file with
a name of the form: "INFILE.pkl". Then, you can substitute this
new file for the corresponding SEGMENTATION or ANNOTATIONS
argument in future Segtools calls and Segtools will parse
it in much faster (but be sure to keep the ".pkl" extension intact)!
"""

from __future__ import division, with_statement

__version__ = "$Revision$"

import sys

from . import Annotations, Segmentation

def preprocess(infile, annotation=False, verbose=False, clobber=False):
    if annotation:
        file_type = Annotations
    else:
        file_type = Segmentation
    print >>sys.stderr, ("Preprocessing into a(n) %s object."
                         " See --help for more options"
                         % file_type.__name__)
    parsed_obj = file_type(infile, verbose=verbose)
    parsed_obj.pickle(verbose=verbose, clobber=clobber)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] INFILE"
    description = __doc__.strip()
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("--clobber", dest="clobber",
                      default=False, action="store_true",
                      help="Overwrite existing output files without warning")
    parser.add_option("-q", "--quiet", action="store_false",
                      dest="verbose", default=True,
                      help="Do not print diagnostic messages.")
    parser.add_option("--annotation", action="store_true",
                      dest="annotation", default=False,
                      help="Read INFILE as a set of annotations, rather"
                      " than as a segmentation (default).")

    (options, args) = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    infile = args[0]
    kwargs = {"verbose": options.verbose,
              "clobber": options.clobber,
              "annotation": options.annotation}
    preprocess(infile, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
