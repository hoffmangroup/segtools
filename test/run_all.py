#!/usr/bin/env python
from __future__ import division, with_statement

import os
import re
import sys
import unittest

def main(verbose=True):
    # Move to test directory to allow imports
    os.chdir(os.path.dirname(__file__) or ".")

    # Gather a list of unittest modules
    filenames = os.listdir(os.getcwd())
    regex = re.compile("^test_.*\.py$", re.IGNORECASE)
    module_filenames = filter(regex.search, filenames)
    make_module_name = lambda filename: filename[:-3]
    modulenames = map(make_module_name, module_filenames)
    if verbose:
        print "Found test modules: %r" % modulenames

    # Run the test suite for each
    suite = unittest.defaultTestLoader.loadTestsFromNames(modulenames)
    if verbose:
        verbosity = 2
    else:
        verbosity = 1

    unittest.TextTestRunner(verbosity=verbosity).run(suite)

if __name__ == "__main__":
    sys.exit(main())
