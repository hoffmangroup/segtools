#!/usr/bin/env python
"""
If Segtools is taking too long to parse your large segmentation, this
tool allows you to pre-process your segmentation once
and have it load faster in the future.
BEDFILE will be parsed to create a special "segmentation" file with
a name of the form: "OUTBASE.<ext>". Then, you can substitute this
new file for the BEDFILE argument in future Segtools calls and Segtools
will parse it in much faster (but be sure to keep the <ext> in tact)!
"""

from __future__ import division, with_statement

__version__ = "$Revision$"

import os
import sys

from . import Segmentation

def pickle_segmentation(infile, outbase, verbose=False, clobber=False):
    segmentation = Segmentation(infile, verbose=verbose)
    segmentation.pickle(outbase, verbose=verbose, clobber=clobber)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] BEDFILE OUTBASE"
    description = __doc__.strip()
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", dest="verbose",
                      default=False, action="store_true",
                      help="Print diagnostic messages")
    parser.add_option("--clobber", dest="clobber",
                      default=False, action="store_true",
                      help="Overwrite existing output files without warning")

    (options, args) = parser.parse_args(args)

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfile = args[0]
    outbase = args[1]
    kwargs = {"verbose": options.verbose,
              "clobber": options.clobber}
    pickle_segmentation(bedfile, outbase, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
