#!/usr/bin/env python

"""segtools: Tools for exploratory analysis of genomic segmentations

Segtools is a Python package designed to put genomic segmentations back
in the context of the genome! Using R for graphics, Segtools provides a
number of modules to analyze a segmentation in various ways and help you
interpret its biological relevance.
"""

__version__ = "1.1.6"

# Copyright 2008-2011 Michael M. Hoffman <mmh1@uw.edu>
# Copyright 2009-2010 Orion J. Buske <stasis@uw.edu>

import sys

# required for from __future__ import division, with_statement;
# relative imports
assert sys.version_info >= (2, 5, 1)

from ez_setup import use_setuptools
use_setuptools()

from setuptools import find_packages, setup

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "http://noble.gs.washington.edu/proj/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
segtools-aggregation = segtools.aggregation:main
segtools-compare = segtools.compare:main
segtools-feature-distance = segtools.feature_distance:main
segtools-flatten = segtools.flatten:main
segtools-gmtk-parameters = segtools.gmtk_parameters:main
segtools-html-report = segtools.html:main
segtools-length-distribution = segtools.length_distribution:main
segtools-nucleotide-frequency = segtools.nucleotide_frequency:main
segtools-overlap = segtools.overlap:main
segtools-preprocess = segtools.preprocess:main
segtools-signal-distribution = segtools.signal_distribution:main
segtools-transition = segtools.transition:main
"""
#segtools = segtools.validate_all:main

install_requires = ["genomedata", "numpy>=1.3", "rpy2", "textinput",
                    "tables>=2.1"]

if __name__ == "__main__":
    setup(name=name,
          version=__version__,
          description=short_description,
          author="Michael Hoffman",
          author_email="mmh1@uw.edu",
          maintainer="Orion Buske",
          maintainer_email="stasis@uw.edu",
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          install_requires=install_requires,
          zip_safe=False,  # For R files to source others, they can't be zip'd
          packages=find_packages("."),
          include_package_data=True,
          entry_points=entry_points
          )
