#!/usr/bin/env python

"""segtools: Tools for exploratory analysis of genomic segmentations

Segtools is a Python package designed to put genomic segmentations back
in the context of the genome! Using R for graphics, Segtools provides a
number of modules to analyze a segmentation in various ways and help you
interpret its biological relevance.
"""

from segtools import __version__

# Copyright 2008-2011 Michael M. Hoffman <mmh1@uw.edu>
# Copyright 2009-2011 Orion J. Buske <stasis@uw.edu>

import sys

# required for from __future__ import division, with_statement;
# relative imports
assert sys.version_info >= (2, 7, 0)

from setuptools import find_packages, setup

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "http://pmgenomics.ca/hoffmanlab/proj/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

classifiers = ["Natural Language :: English",
               "Development Status :: 5 - Production/Stable",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GNU General Public License v2 "
               "(GPLv2)",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering :: Bio-Informatics",
               "Operating System :: Unix",
               "Programming Language :: Python :: 2.7",
               "Programming Language :: Python :: 3",
               "Programming Language :: Python :: 3.6"
               ]


def entry_points():
    python3_compatible_entry_points = [
        "segtools-aggregation = segtools.aggregation:main",
        "segtools-feature-distance = segtools.feature_distance:main",
        'segtools-gmtk-parameters = segtools.gmtk_parameters:main',
        'segtools-length-distribution = segtools.length_distribution:main',
        "segtools-preprocess = segtools.preprocess:main"
    ]

    python2_only_entry_points = [
        "segtools-compare = segtools.compare:main",
        "segtools-flatten = segtools.flatten:main",
        "segtools-html-report = segtools.html:main",
        "segtools-nucleotide-frequency = segtools.nucleotide_frequency:main [genomedata]",
        "segtools-overlap = segtools.overlap:main",
        "segtools-relabel = segtools.relabel:main",
        "segtools-signal-distribution = segtools.signal_distribution:main [genomedata]",
        "segtools-transition = segtools.transition:main"
    ]

    entry_points = {
            'console_scripts': python3_compatible_entry_points
        }

    if sys.version_info[0] == 2:
        entry_points['console_scripts'].extend(python2_only_entry_points)

    return entry_points


install_requires = ["numpy>=1.3", "rpy2>=2.6.0,<2.9", "six", "textinput"]
# XXX: add optional requirement for PyGraphviz
extras_require = {'genomedata': "genomedata"}

if __name__ == "__main__":
    setup(name=name,
          version=__version__,
          description=short_description,
          author="Michael Hoffman",
          author_email="michael.hoffman@utoronto.ca",
          maintainer="Michael Hoffman",
          maintainer_email="michael.hoffman@utoronto.ca",
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          install_requires=install_requires,
          extras_require=extras_require,
          zip_safe=False,  # For R files to source others, they can't be zip'd
          packages=find_packages("."),
          package_data={name: ['R/*.R',
                               'resources/*.tmpl']},
          entry_points=entry_points()
          )
