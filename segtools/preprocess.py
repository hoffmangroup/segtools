#!/usr/bin/env python
"""
If Segtools is taking too long to parse your large segmentation or
annotation file (specify --annotation), this
tool allows you to pre-process that file only once
and have it load much faster in the future.
INFILE will be parsed to create a special binary file with
a name of the form: "INFILE.pkl" or "OUTFILE.pkl".
Then, you can substitute this
new file for the corresponding SEGMENTATION or ANNOTATION
argument in future Segtools calls and Segtools will parse
it in much faster (but be sure to keep the ".pkl" extension intact)!
"""

from __future__ import division, with_statement

import sys

from . import Annotation, Segmentation, add_common_options, log
from .version import __version__

def preprocess(infilename, outfilename=None, annotation=False, verbose=False, clobber=False):
    if annotation:
        file_type = Annotation
    else:
        file_type = Segmentation

    log("Making %s object. See --help for more options" % file_type.__name__,
        verbose=verbose)
    parsed_obj = file_type(infilename, verbose=verbose)
    parsed_obj.pickle(outfilename, verbose=verbose, clobber=clobber)

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTIONS] INFILE [OUTFILE]"
    description = __doc__.strip()
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    add_common_options(parser, ['clobber', 'quiet'])
    parser.add_option("--annotation", action="store_true",
                      dest="annotation", default=False,
                      help="Read INFILE as an annotation, rather"
                      " than as a segmentation (default).")

    (options, args) = parser.parse_args(args)

    if len(args) < 1 or len(args) > 2:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    kwargs = dict(options.__dict__)
    preprocess(*args, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
