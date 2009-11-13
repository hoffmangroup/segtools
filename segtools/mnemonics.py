#!/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
mnemonics.py

"""

import os
import sys

from rpy2.robjects import r, rinterface, numpy2ri

from .common import check_clobber, make_tabfilename, r_source

def create_mnemonic_file(datafile, dirpath, clobber=False,
                         namebase=None, gmtk=False):
    """Generate a mnemonic file with R clustering code.

    Datafile can either be a signal stats file or a GMTK params file
      (in which case gmtk=True must be specified to make.mnemonic.file)
    Calls R code that writes mnemonics to a file.
    Returns name of created mnemonic file
    """
    assert os.path.isfile(datafile)

    r_source("common.R")
    r_source("track_statistics.R")

    if namebase is None:
        namebase = os.path.basename(datafile)
        if namebase.endswith(".tab"):
            namebase = namebase[:-4]
        if namebase.endswith(".params"):
            namebase = namebase[:-7]

    mnemonic_base = os.extsep.join([namebase, "mnemonics"])

    filename = make_tabfilename(dirpath, mnemonic_base)
    check_clobber(filename, clobber)

    # Create mnemonic file via R
    try:
        r_filename = r["make.mnemonic.file"](datafile, filename,
                                             gmtk=gmtk)
        mnemonicfilename = str(r_filename[0])

        print >>sys.stderr, "Generated mnemonic file: %s" % mnemonicfilename
        return mnemonicfilename
    except rinterface.RRuntimeError:
        print >>sys.stderr, ("ERROR: Failed to create mnemonic file."
                             " Continuing without mnemonics.")
        return None

if __name__ == "__main__":
    pass
