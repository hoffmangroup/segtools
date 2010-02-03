================================
Segtools |version| documentation
================================
:Homepage: http://noble.gs.washington.edu/proj/segtools
:Author: Orion Buske <stasis at uw dot edu>
:Organization: University of Washington
:Address: Department of Genome Sciences, PO Box 355065, 
          Seattle, WA 98195-5065, United States of America
:Copyright: 2009, Orion Buske
:Last updated: |today| 

.. currentmodule:: segtools

Description
===========
Segtools is a python package designed to put genomic segmentations back
in the context of the genome! Using R for graphics, segtools provides a
number of modules to analyze a segmentation in various ways and help
you interpret its biological relevance.

Segmentations should be in ``BED`` format, with the ``name`` field of each
line used specifying the segment label of that line. The segtools modules 
allow you to compare the properties of the segment labels with one another.

Don't hesistate to contact the author for more information!


Installation
============
A simple, interactive script_ has been created to install segtools 
(and most dependencies) on any Linux platform. Installation is as simple
as downloading and running this script! For instance::
   
   wget http://noble.gs.washington.edu/proj/segtools/install.py
   python install.py

.. _script: http://noble.gs.washington.edu/proj/segtools/install.py

.. note:: 
   The following are prerequisites:
          
   - Python 2.5.1+
   - Zlib


Usage
=====

All segtools modules requires, at the very least, a segmentation in 
`BED format`_. Some modules also require additional files, such as 
genomic feature files to compare with the segmentation. Modules are most
easily used from the command line interface, but all can be loaded and 
run straight from python if desired (though this is not yet documented).

.. _`BED format`: http://genome.ucsc.edu/FAQ/FAQformat#format1

The following commands are currently available:
    
:program:`segtools-feature-aggregation`
      Analyzes the relative occurrance of each segment label 
      around the provided genomic features.
:program:`segtools-label-transition`
      Analyses the transitions between segment labels and
      the structure of their interaction.
:program:`segtools-length-distribution`
      Analyzes the distribution of segment lengths
      and their coverage of the genome for each segment label.
:program:`segtools-signal-distribution`
      Analyzes the distribution of genomic signal 
      tracks for each segment label.
:program:`segtools-nucleotide-frequency`
      Analyzes the frequencies of nucleotides and
      dinucleotides in segments of each label.
:program:`segtools-overlap`
      Analyzes the frequency with which each segment label overlaps
      features of various types.
:program:`segtools-html-report`
      Combines the output of the other modules and generates
      an html report for viewing.

All the above commands respond to ``-h`` and ``--help``.

Each module generates:

     - tab-delimited (``tab``) data files
     - image files (in ``png`` and ``pdf`` format and in 
       normal, thumbnail, and slide layouts), and 
     - partial HTML (``div``) files.



.. Technical description
.. ---------------------


Modules
=======


:mod:`.feature_aggregation`
---------------------------

.. module:: segtools.feature_aggregation

This modules aggregates segmentation data around features, generating
a histogram for each segmentation label that shows the frequency of
observing that label at that position relative the the feature.

If using gene mode, the input file should have features with names:
exon, start_codon, CDS
as provided by exporting data from the `UCSC Table Browser`_ in GFF format.

.. _`UCSC Table Browser`: http://genome.ucsc.edu/cgi-bin/hgTables?command=start




:mod:`.html`
-------------------

.. module:: segtools.html

This module is intended to be run after other segtools modules. It searches
the local (or provided) directory for ``div`` files produced by the
other segtools modules and compiles the data into an HTML report for 
review.

.. program:: segtools-html-report


The ``BEDFILE`` argument and :option:`--mnemonic-file` option 
should be the same as used to run the other segtools modules.


:mod:`.label_transition`
------------------------

.. module:: segtools.label_transition

Provides command-line and package entry points for analyzing the observed
segmentation label transitions in the given BED-formatted segmentation.


Accepts an input file containing a matrix of transition
probabilities and generates several output files:

  - a heatmap of the matrix
  - a graph of a thresholded form of the transition matrix

.. program:: segtools-label-transition



:mod:`.length_distribution`
---------------------------

.. module:: segtools.length_distribution

Provides command-line and package entry points for analyzing the segment
length distribution in a provided BED-formatted segmentation.

.. program:: segtools-length-distribution



:mod:`.nucleotide_frequency`
----------------------------

.. module:: segtools.nucleotide_frequency

Provides command-line and package entry points for analyzing nucleotide
and dinucleotide frequencies for each segmentation label.

.. program:: segtools-nucleotide-frequency



:mod:`.overlap`
---------------

.. module:: segtools.overlap

Evaluates the overlap between two BED files, based upon a specification 
that can be found here__.

__ http://encodewiki.ucsc.edu/EncodeDCC/index.php/
   Overlap_analysis_tool_specification

.. program:: segtools-overlap



:mod:`.signal_distribution`
---------------------------

.. module:: segtools.signal_distribution

Provides command-line and package entry points for analyzing the signal
distribution over tracks and labels.

.. program:: segtools-signal-distribution
