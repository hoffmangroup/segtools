#!/usr/bin/env python

__version__ = "$Revision: 84 $"

# Copyright 2009 Orion Buske <orion.buske@gmail.com>

import os
import re
import sys

from pkg_resources import resource_listdir

def test(verbose=False):
    import unittest

    # Gather a list of unittest modules
    filenames = resource_listdir(__name__, ".")
    regex = re.compile("^test_.*\.py$", re.IGNORECASE)
    module_filenames = filter(regex.search, filenames)
    def make_module_name(filename):
        return os.extsep.join([__name__, filename[:-3]])
    modulenames = map(make_module_name, module_filenames)
    print "Found test modules: %r" % modulenames
    map(__import__, modulenames)  # Import them all
    modules = [sys.modules[modulename] for modulename in modulenames]
    # Run the test suite for each
    suite = unittest.TestSuite([module.suite() for module in modules])
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    unittest.TextTestRunner(verbosity=verbosity).run(suite)

if __name__ == "__main__":
    pass
