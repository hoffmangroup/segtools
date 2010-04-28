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
   For this release (|release|), the documentation is *still* very
   incomplete. This will be remedied as soon as possible. In the mean
   time, please don't hesitate to contact the author (above) with
   *ANY* questions you might have concerning the installation, use,
   and interpretation of results of the Segtools package. We think
   this tool is very useful for a wide array of applications, but the user
   interface and documentation are a little rough at the moment. Thank
   you for your patience as we try to refine this package.


Installation
============
A simple, interactive script_ has been created to install Segtools 
(and most dependencies) on any Linux platform. Installation is as simple
as downloading and running this script! For instance::
   
   wget http://noble.gs.washington.edu/proj/segtools/install.py
   python install.py

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

If you want to change the order in which labels appear or the
text displayed in plots, a `mnemonic file`_ can be created.
Segtools commands can the be re-run with the ``--replot`` flag and
the ``--mnemonic-file=<FILE>`` option to regenerate the plots without
redoing the computation. Similarly, mnemonic files can be swapped or
revised and new images created with relative ease.

Most Segtools commands look for patterns between segment labels in a
segmentation and some known annotation. For such commands, the annotations
should usually specified in `GFF format` (although some commands accept
BED files as well and some require GTF_ or Genomedata_ formats).

.. Workflow
.. ========
.. Coming soon.

Usage
=====

All Segtools commands require, at the very least, a segmentation in 
`BED format`_. Some commands also require additional files, such as 
genomic feature files to compare with the segmentation. It is easiest
to run Segtools through its command line interface, but the same
functionality can be achieved through its Python interface
(though this is not yet documented).

Command-line interface
----------------------

Core commands:

- :ref:`segtools-feature-aggregation <segtools-feature-aggregation>`: 
  Analyzes the relative 
  occurrance of each segment label around the provided genomic features.
- :ref:`segtools-label-transition <segtools-label-transition>`: 
  Analyzes the transitions 
  between segment labels and the structure of their interaction.
- :ref:`segtools-length-distribution <segtools-length-distribution>`: 
  Analyzes the distribution 
  of segment lengths and their coverage of the genome for each segment label.
- :ref:`segtools-signal-distribution <segtools-signal-distribution>`:
  Analyzes the distribution 
  of genomic signal tracks for each segment label.
- :ref:`segtools-nucleotide-frequency <segtools-nucleotide-frequency>`: 
  Analyzes the frequencies 
  of nucleotides and dinucleotides in segments of each label.
- :ref:`segtools-overlap <segtools-overlap>`: 
  Analyzes the frequency with which 
  each segment label overlaps features of various types.
- :ref:`segtools-html-report <segtools-html-report>`: 
  Combines the output of the 
  other commands and generates an html report for easier viewing.

Utility commands:

- :ref:`segtools-bed-compare <segtools-bed-compare>`: 
  Measure base-wise edit distance 
  between two segmentations.
- :ref:`segtools-feature-distance <segtools-feature-distance>`: 
  Reports the distance from 
  each segment to the nearest feature in each of a list of feature files.
- :ref:`segtools-flatten-bed <segtools-flatten-bed>`: 
  General tool for flattening 
  overlapping segments, but flattens them into segments defined 
  by the set of segment labels that overlap the region. 

Other commands:

- :ref:`segtools-gmtk-parameters <segtools-gmtk-parameters>`: 
  Analyzes GMTK_ emission parameters and state transitions.


All the above commands respond to ``-h`` and ``--help``, and this will
display the most up-to-date usage information and options.

Where relevant, commands accept `mnemonic files <mnemonic file>`_
through the ``--mnemonic-file`` option.

Each command generates:

     - tab-delimited (``tab``) data files
     - image files (in ``png`` and ``pdf`` format and in 
       normal, thumbnail, and slide layouts), and 
     - partial HTML (``div``) files.

.. Technical description
.. ---------------------


Commands
========


.. _segtools-bed-compare:

segtools-bed-compare
--------------------

.. program:: segtools-bed-compare




.. ####################### FEATURE AGGREGATION #######################

.. _segtools-feature-aggregation:

segtools-feature-aggregation
----------------------------

.. program:: segtools-feature-aggregation

This command looks at the aggregate occurance of segment labels around
and within annotated features. A typical example of this would be to
look at the relative occurances of segment labels around transcription
start sites (TSSs). You would do this with something like::

  segtools-feature-aggregation --normalize segmentation.bed tss.gff

If you had two different classes of TSSs that you were interrested in
(say, expressed and unexpressed), you can use the 3rd column of the GFF
file as a grouping variable and then specify the :option:`--groups` flag. 

By default, the y-axis of the aggregation plot is the number of 
segments of a given label that overlap a region. This is useful 
in some applications, but more often you are interested in the
enrichment or depletion of various labels relative to what you 
might expect by chance. This is especially true if the segments in
one label are significantly longer than those in another label.
In these cases, the :option:`--normalize` flag should be used.

**Significance**:
  Whether the y-axis is normalized or not, the significance of the
  overlap at a region is shown in the plot. If :option:`--groups` is not
  specified or there is only one group, then significance is shown by
  shading the regions that are significant. Otherwise, the significance
  of the various groups are shown using colored "rugs" at the bottom of 
  the plot. The probability of observing :math:`n` overlapping segments
  of a given label at a given position is modeled with a binomial
  distribution: :math:`p = binom(n, N, f_{rand})`, where :math:`N` is the
  total number of overlapping segments at that position and
  :math:`f_{rand}` is the same as in :option:`--normalize`.
  The p-value is then the probability of observing an overlap count
  as extreme or more extreme than :math:`n` (either enrichment or
  depletion). This corresponds to a two-tailed binomial test.

  .. Didn't use note directive because of latex math image backgrounds.

  *Note:* If :math:`n > 100` and :math:`f_{rand}*N < 10`, a Poisson
  approximation is used.


**Selected options**:

.. cmdoption:: --help, -h

   Display complete usage information

.. cmdoption:: --mode <mode>, -m <mode>

   Specify the aggregation mode. The following options are available: 
   ``point``, ``region``, and ``gene``. The default is ``point``.

   ``point``: 
   This mode aggregates around point-like features such as TSSs, TESs,
   or single-base peak calls. This mode looks at where segments 
   occur in the 5' and 3' flanking regions of each feature. If the
   feature annotations have strand specifications (column 7), the aggregation
   is strand-corrected so that the 5' flank region is always upstream
   of the feature. The width (in base pairs) of these flanking regions 
   can be set with the ``--flank-bases`` option (default 500 bp).

   ``region``:
   This mode aggregates around region-like features such as
   transcription factor binding sites, ChIP-seq peak calls,
   or promoter regions. This will be the appropriate mode for most
   annotations. This mode is similar to ``point``, but with the addition
   of an ``internal`` region which is aggregated over as well. To account
   for regions of varying length, evenly-spaced samples are taken from
   the span of each feature. The number of these samples can be set
   with ``--region-samples``. Features than span fewer bases than this
   sample number are skipped.

   ``gene``:
   This is a special mode for aggregating with respect to an idealized
   gene model. Rather than specifying a normal GFF file, the annotation file
   must be in `GTF format`_ and have features with names: ``exon``, ``CDS`,`
   as provided by exporting data from the `UCSC Table Browser`_ in 
   `GTF format`_. This mode is similar to ``region``, but with many regions
   that correspond to idealized transcriptional and translational models of
   genes. For the transcriptional model, there are regions
   corresponding to initial, internal, and terminal exons and introns. 
   For the translational model, there are initial, internal, and terminal
   5' and 3' UTR regions, and initial and terminal CDSs. These two models
   are layed out in logical progressions so that genes are
   viewed in profile and gene-component-specific associations can
   be easily seen. Because introns and exons
   are typically different lengths, ``--intron-samples`` and
   ``--exon-samples`` options allow you to specify the number of
   samples that are taken in these regions (like in ``region`` mode).
   *Note: If there are multiple transcripts with the same
   gene ID, the longest transcript is used.*

.. cmdoption:: --normalize

   This option normalizes the y-axis of the aggregation plot,
   displaying enrichment and depletion instead of counts at each
   position. The enrichment of label :math:`l` at position :math:`p`
   is calulated with the following formula:

   .. math:: enrichment(l, p) = \log_2 \dfrac{f_{obs} + 1}{f_{rand} + 1}

   where :math:`f_{obs}` is the observed overlap frequency and
   :math:`f_{obs}` is the overlap frequency expected at random, 
   defined by:
   
   .. math:: 

      f_{obs} = \dfrac{count(l, p)}{\sum_{labels} count(p)}

      f_{rand} = \dfrac{bases\_in\_label(l)}{\sum_{labels} bases\_in\_label}

   The enrichment is thus bounded by :math:`[-1,1]`, with 1
   corresponding to extreme enrichment and -1 to extreme depletion.

.. cmdoption:: --groups

   Group the features by the value of the 3rd column of the GFF or GTF
   file (the ``name`` field). This is useful if you wanted to compare
   the aggregation profiles with respect to multiple classes of
   features, such as TSSs split by expression level or cell type.





.. ####################### FEATURE DISTANCE #######################

.. _segtools-feature-distance:

segtools-feature-distance
-------------------------

.. program:: segtools-feature-distance





.. ####################### FLATTEN BED #######################

.. _segtools-flatten-bed:

segtools-flatten-bed
--------------------

.. program:: segtools-flatten-bed




.. ####################### GMTK PARAMETERS #######################

.. _segtools-gmtk-parameters:

segtools-gmtk-parameters
------------------------

.. program:: segtools-gmtk-parameters




.. ####################### HTML REPORT #######################

.. _segtools-html-report:

segtools-html-report
--------------------

This command is intended to be run after other Segtools commands. It searches
the local (or provided) directory for ``div`` files produced by the
other Segtools commands and compiles the data into an HTML report for 
review.

.. program:: segtools-html-report


The ``BEDFILE`` argument and ``--mnemonic-file`` option 
should be the same as used to run the other Segtools commands.



.. ####################### LABEL TRANSITION #######################

.. _segtools-label-transition:

segtools-label-transition
-------------------------

.. program:: segtools-label-transition




.. ####################### LENGTH DISTRIBUTION #######################

.. _segtools-length-distribution:

segtools-length-distribution
----------------------------

.. program:: segtools-length-distribution



.. ####################### NUCLEOTIDE FREQUENCY #######################

.. _segtools-nucleotide-frequency:

segtools-nucleotide-frequency
-----------------------------

.. program:: segtools-nucleotide-frequency





.. ####################### OVERLAP #######################

.. _segtools-overlap:

segtools-overlap
----------------

.. program:: segtools-overlap




.. ####################### SIGNAL DISTRIBUTION #######################

.. _segtools-signal-distribution:

segtools-signal-distribution
----------------------------

.. program:: segtools-signal-distribution






.. _`mnemonic file`:

Mnemonics
=========

Mnemonic files are supported by most of the Segtools commands and provide
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
into Segtools commands with ``--mnemonic-file=/path/to/second_try.mnemonics``.

If you had previously run Segtools commands on the segmentation before
creating these mnemonics, you could speed up the plot corrections by 
using the command's ``--replot`` option (all other options and
arguments should still be specified to ensure correctness).


.. _`BED4+ format`:
.. _`BED format`: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _script: http://noble.gs.washington.edu/proj/segtools/install.py

.. _`GFF format`: http://genome.ucsc.edu/FAQ/FAQformat.html#format3

.. _GTF:
.. _`GTF format`: http://genome.ucsc.edu/FAQ/FAQformat.html#format4
.. _GMTK: http://ssli.ee.washington.edu/~bilmes/gmtk/
.. _Genomedata: http://noble.gs.washington.edu/proj/genomedata/
.. _`UCSC Table Browser`: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
