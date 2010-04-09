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

..  Buske OJ, Hoffman MM, Noble WS, "Exploratory analysis of genomic
..  segmentations with Segtools." In preparation.

Segtools is a Python package designed to put genomic segmentations back
in the context of the genome! Using R for graphics, Segtools provides a
number of modules to analyze a segmentation in various ways and help
you interpret its biological relevance. These modules have command-line
and Python interfaces, but the command-line interfaces are most thoroughly
documented.

.. warning::
   For this release (1.0.0), the documentation is *very* incomplete. This
   will be remedied as soon as possible. In the mean time, please don't
   hesitate to contact the author (above) with *ANY* questions 
   you might have concerning the installation, use, and interpretation of 
   results of the Segtools package. We think this tool is very 
   useful for a wide array of applications, but the user
   interface and documentation are a little rough at the moment. Thank
   you for your patience as we try to refine this package.


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

Basics
======
A segmentation is typically a partition of a genome (or part of a genome)
into non-overlapping segments, each of which is assigned one of a small
set of labels. The idea is that segments that share a common label are somehow
similar, and those that have different labels are somehow different.
Segtools helps you identify the similarities and differences between these
labels to help you understand your segmentation at a higher level.

Input
=====
Segmentations should be in `BED4+ format`_, with one 
segment per line and the ``name`` field used to specify the segment label. 
For best results, the number of unique segment labels should be between
2 and around 40.

.. _`BED4+ format`: http://genome.ucsc.edu/FAQ/FAQformat#format1

If you want to change the order in which labels appear or the
text displayed in plots, a `mnemonic file`_ can be created.
Segtools commands can the be re-run with the ``--replot`` flag and
the ``--mnemonic-file=<FILE>`` option to regenerate the plots without
redoing the computation. Similarly, mnemonic files can be swapped or
revised and new images created with relative ease.

Most Segtools modules look for patterns between segment labels in a
segmentation and some known annotation. For such modules, the annotations
should usually specified in `GFF format` (although some modules accept
BED files as well and some require GTF_ or Genomedata_ formats).

.. _`GFF format`: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
.. _GTF: http://genome.ucsc.edu/FAQ/FAQformat.html#format4
.. _Genomedata: http://noble.gs.washington.edu/proj/genomedata/

Workflow
========
Coming soon.

Usage
=====

All Segtools modules require, at the very least, a segmentation in 
`BED format`_. Some modules also require additional files, such as 
genomic feature files to compare with the segmentation. Modules are most
easily run from the command line, but all can be loaded and 
run straight from python if desired (though this is not yet documented).

.. _`BED format`: http://genome.ucsc.edu/FAQ/FAQformat#format1

Command-line interface
----------------------

Core commands:

- :program:`segtools-feature-aggregation`: Analyzes the relative 
  occurrance of each segment label around the provided genomic features.
- :program:`segtools-label-transition`: Analyzes the transitions 
  between segment labels and the structure of their interaction.
- :program:`segtools-length-distribution`: Analyzes the distribution 
  of segment lengths and their coverage of the genome for each segment label.
- :program:`segtools-signal-distribution`: Analyzes the distribution 
  of genomic signal tracks for each segment label.
- :program:`segtools-nucleotide-frequency`: Analyzes the frequencies 
  of nucleotides and dinucleotides in segments of each label.
- :program:`segtools-overlap`: Analyzes the frequency with which 
  each segment label overlaps features of various types.
- :program:`segtools-html-report`: Combines the output of the 
  other modules and generates an html report for easier viewing.

Utility commands:

- :program:`segtools-bed-compare`: Measure base-wise edit distance 
  between two segmentations.
- :program:`segtools-feature-distance`: Reports the distance from 
  each segment to the nearest feature in each of a list of feature files.
- :program:`segtools-flatten-bed`: General tool for flattening 
  overlapping segments, but flattens them into segments defined 
  by the set of segment labels that overlap the region. 

Other commands:

- :program:`segtools-gmtk-parameters`: Analyzes GMTK_ emission 
  parameters and interfaces with appropriate Segtools modules.

.. _GMTK: http://ssli.ee.washington.edu/~bilmes/gmtk/

All the above commands respond to ``-h`` and ``--help``.

Each command generates:

     - tab-delimited (``tab``) data files
     - image files (in ``png`` and ``pdf`` format and in 
       normal, thumbnail, and slide layouts), and 
     - partial HTML (``div``) files.

.. Technical description
.. ---------------------


Modules
=======

.. _segtools-feature-aggregation:

:mod:`.feature_aggregation`
---------------------------

This modules aggregates segment labels around features, generating
a histogram for each segmentation label that shows the likelihood of
observing that label at that position relative to the feature.

If using ``gene`` mode, the input file should have features with names:
``exon``, ``start_codon``, ``CDS``
as provided by exporting data from the `UCSC Table Browser`_ in `GFF format`_.

.. _`UCSC Table Browser`: http://genome.ucsc.edu/cgi-bin/hgTables?command=start


.. program:: segtools-feature-aggregation


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


.. program:: segtools-label-transition



:mod:`.length_distribution`
---------------------------

.. module:: segtools.length_distribution


.. program:: segtools-length-distribution



:mod:`.nucleotide_frequency`
----------------------------

.. module:: segtools.nucleotide_frequency



.. program:: segtools-nucleotide-frequency



:mod:`.overlap`
---------------

.. module:: segtools.overlap


.. program:: segtools-overlap



:mod:`.signal_distribution`
---------------------------

.. module:: segtools.signal_distribution


.. program:: segtools-signal-distribution


.. _`mnemonic file`:

Mnemonic file
=============
Mnemonic files are supported by all the Segtools modules and provide
a way to rename and reorder the displayed labels without repeating
the full analysis. Mnemonic files must be two or three tab-separated 
columns and must contain start with the following header 
(the description column is optional)::

  old{TAB}new[{TAB}description]
  
**Renaming**:

Each line of the mnemonic file specifies a mapping from the "old"
name (the one appearing in the segmentation file) to the "new" name
(the one to be displayed). Since the new name must fit into small spaces
in plots (such as axis labels), it is recommended for this field to be
a few characters (such as "I" for insulator). Longer descriptions can
be specified in the description column.

**Reordering**:

The order of the lines in the mnemonic file determines the order
the labels will be shown in plots.

**Example**:

If the segmentation file contains segments with labels of ``A``, ``B``,
and ``C``, but realized you wanted ``A`` to be displayed as ``A1``,
``C`` to be displayed as ``A2``, and the two of them to be next to
each other in plots, you should construct the following mnemonic file::

  old	new
  A	A1
  C	A2
  B	B

Including the B line is not neccessary, but it makes it easier to
reorder the labels later (for instance, if you want B to come first).
A description column could also have been included. This file should be
saved as something like ``second_try.mnemonics`` and should be passed
into Segtools modules with ``--mnemonic-file=/path/to/second_try.mnemonics``.

If you had previously run Segtools modules on the segmentation before
creating these mnemonics, you could speed up the plot corrections by 
using the module's ``--replot`` option.

.. warning:: 
   Mnemonics are almost completely supported, with the exception of the
   tables appearing in the HTML report. The names in these tables will
   correspond to the names used during the first run of the module
   and are not updated with ``--replot``. For those modules that generate
   a table, the full analysis must be redone to regenerate these tables.
   We hope to fix this soon.