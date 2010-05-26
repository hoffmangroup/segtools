#!/usr/bin/env python
"""
If Segtools is taking too long to parse your large segmentation or
annotation file, this
tool allows you to pre-process your segmentation once
and have it load faster in the future.
INFILE will be parsed to create a special binary file with
a name of the form: "INFILE.pkl". Then, you can substitute this
new file for the SEGMENTATION or ANNOTATIONS argument in future
Segtools calls and Segtools will parse it in much faster
(but be sure to keep the ".pkl" extension intact)!
"""

from __future__ import division, with_statement

__version__ = "$Revision$"

import sys

from . import Segmentation

def preprocess(infile, verbose=False, clobber=False):
    segmentation = Segmentation(infile, verbose=verbose)
    segmentation.pickle(verbose=verbose, clobber=clobber)

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

    (options, args) = parser.parse_args(args)

    if len(args) != 1:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    infile = args[0]
    kwargs = {"verbose": options.verbose,
              "clobber": options.clobber}
    preprocess(infile, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
