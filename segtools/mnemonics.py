#!/bin/env python
from __future__ import division, with_statement

"""
mnemonics.py

"""

import os

from . import log, RInterface
from .common import check_clobber

def create_mnemonic_file(datafile, dirpath, clobber=False,
                         namebase=None, gmtk=False, verbose=True):
    """Generate a mnemonic file with R clustering code.

    Datafile can either be a signal stats file or a GMTK params file
      (in which case gmtk=True must be specified to make.mnemonic.file)
    Calls R code that writes mnemonics to a file.
    Returns name of created mnemonic file
    """
    assert os.path.isfile(datafile)

    R = RInterface(["common.R", "track_statistics.R"], verbose=verbose)

    if namebase is None:
        namebase = os.path.basename(datafile)
        if namebase.endswith(".tab"):
            namebase = namebase[:-4]
        if namebase.endswith(".params"):
            namebase = namebase[:-7]

    mnemonic_base = os.extsep.join([namebase, "mnemonics"])

    filename = os.path.join(dirpath, mnemonic_base)
    check_clobber(filename, clobber)

    # Create mnemonic file via R
    try:
        r_filename = R.call("make.mnemonic.file",
                            datafile, filename, gmtk=gmtk)
        # Peel off rpy2 layers
        mnemonicfilename = str(r_filename[0])

        log("Created mnemonic file: %s" % mnemonicfilename, verbose)

        return mnemonicfilename
    except R.RError:
        log("ERROR: Failed to create mnemonic file."
            " Continuing without mnemonics.")
        return None

if __name__ == "__main__":
    pass
