================================
Segtools |version| documentation
================================
:Author: Orion Buske <stasis at uw dot edu>
:Organization: University of Washington
:Address: Department of Genome Sciences, PO Box 355065, 
          Seattle, WA 98195-5065, United States of America
:Copyright: 2009 Orion Buske
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


Command-line usage summary
..........................

::

 Usage: segtools-feature-aggregation [OPTIONS] BEDFILE FEATUREFILE

 FEATUREFILE should be in GFF or GTF format

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit

   Input options:
     --mnemonic-file=MNEMONICFILENAME
                         If specified, labels will be shown using mnemonics
                         found in this file
     -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: feature_aggregation]

   Flags:
     --clobber           Overwrite existing output files if the specified
                         directory already exists.
     --quick             Compute values only for one chromosome.
     --replot            Load data from output tab files and regenerate plots
                         instead of recomputing data
     --noplot            Do not generate plots
     --no-factor         Do not separate data into different factors

   Aggregation options:
     -m MODE, --mode=MODE
                         one of: ['point', 'region', 'gene'], --gene not
                         implemented [default: point]
     -f FLANKBINS, --flank-bins=FLANKBINS
                         Aggregate this many base pairs off each end of feature
                         or gene [default: 500]
     -r REGIONBINS, --region-bins=REGIONBINS
                         If --mode=region, aggregate over each internalfeature
                         using this many evenly-spaced bins [default: 50]
     -i INTRONBINS, --intron-bins=INTRONBINS
                         If --mode=gene, Aggregate over each intronusing this
                         many evenly-spaced bins [default: 50]
     -e EXONBINS, --exon-bins=EXONBINS
                         If --mode=gene, Aggregate over each exonusing this
                         many evenly-spaced bins [default: 25]


:mod:`.html`
-------------------

.. module:: segtools.html

This module is intended to be run after other segtools modules. It searches
the local (or provided) directory for ``div`` files produced by the
other segtools modules and compiles the data into an HTML report for 
review.

.. program:: segtools-html-report

Command-line usage summary
..........................

::

 Usage: segtools-html-report [OPTIONS] BEDFILE

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   --clobber             Overwrite existing output files if the specified
                         directory already exists.
   --mnemonic-file=MNEMONICFILE
                         If specified, this mnemonic mapping will be included
                         in the report (this should be the same mnemonic file
                         used by the individual modules)
   --results-dir=RESULTSDIR
                         This should be the directory containing all the module
                         output directories (`ls` should return things like
                         "length_distribution/", etc) [default: .]
   -o OUTFILE, --outfile=OUTFILE
                         HTML report file (must be in current directory
                         [default: index.html]

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

Command-line usage summary
..........................

::

 Usage: segtools-label-transition [OPTIONS] BEDFILE

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   --clobber             Overwrite existing output files if the specified
                         directory already exists.
   --noplot              Do not generate transition plots
   --nograph             Do not generate transition graph
   --mnemonic-file=MNEMONIC_FILE
                         If specified, labels will be shown using mnemonics
                         found in this file
   -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: label_transition]
   --gmtk-params=GMTK_FILE
                         If specified, parameters in the given GMTK file will
                         be used to generate plots instead of the observed
                         transitions in the BEDFILE. The BEDFILE will not be
                         used

   Transition frequency plot options:
     --dd, --dendrogram  include dendrogram along edge of levelplot [default:
                         False]

   Transition graph options:
     -p P_THRESH, --prob-threshold=P_THRESH
                         ignore all transitions with probabilities below this
                         absolute threshold [default: 0.15]
     -q Q_THRESH, --quantile-threshold=Q_THRESH
                         ignore transitions with probabilities below this
                         probability quantile [default: 0.0]

   Non-segmentation files:
     --gmtk-params=GMTK_FILE
                         If specified, parameters in the given GMTK file will
                         be used to generate plots instead of the observed 
                         transitions in the BEDFILE. The BEDFILE will not be 
                         used

:mod:`.length_distribution`
---------------------------

.. module:: segtools.length_distribution

Provides command-line and package entry points for analyzing the segment
length distribution in a provided BED-formatted segmentation.

.. program:: segtools-length-distribution

Command-line usage summary
..........................

::

 Usage: segtools-length-distribution [OPTIONS] BEDFILE

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   --clobber             Overwrite existing output files if the specified
                         directory already exists.
   --replot              Load data from output tab files and regenerate plots
                         instead of recomputing data
   --noplot              Do not generate plots
   --mnemonic-file=MNEMONICFILENAME
                         If specified, labels will be shown using mnemonics
                         found in this file
   -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: length_distribution]




:mod:`.nucleotide_frequency`
----------------------------

.. module:: segtools.nucleotide_frequency

Provides command-line and package entry points for analyzing nucleotide
and dinucleotide frequencies for each segmentation label.

.. program:: segtools-nucleotide-frequency

Command-line usage summary
..........................

::

 Usage: segtools-nucleotide-frequency [OPTIONS] BEDFILE GENOMEDATADIR

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   --clobber             Overwrite existing output files if the specified
                         directory already exists.
   --quick               Compute values only for one chromosome.
   --replot              Load data from output tab files and regenerate plots
                         instead of recomputing data
   --noplot              Do not generate plots
   --mnemonic-file=MNEMONICFILENAME
                         If specified, labels will be shown using mnemonics
                         found in this file
   -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: nucleotide_frequency]


:mod:`.overlap`
---------------

.. module:: segtools.overlap

Evaluates the overlap between two BED files, based upon a specification 
that can be found here__.

__ http://encodewiki.ucsc.edu/EncodeDCC/index.php/
   Overlap_analysis_tool_specification

.. program:: segtools-overlap

Command-line usage summary
..........................

::

 Usage: segtools-overlap [OPTIONS] BEDFILE FEATUREFILE

 BEDFILE and FEATUREFILE should both be in BED3+ format (gzip'd okay). BEDFILE
 should correspond to a segmentation. Overlap analysis will be performed in
 both directions (BEDFILE as SUBJECTFILE and QUERYFILE). See for full
 specification: http://encodewiki.ucsc.edu/EncodeDCC/index.php/Overlap_analysis
 _tool_specification

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit

   Flags:
     --clobber           Overwrite existing output files if the specified
                         directory already exists.
     --quick             Compute values only for one chromosome.
     --replot            Load data from output tab files and regenerate plots
                         instead of recomputing data
     --noplot            Do not generate plots

   Parameters:
     -b BY, --by=BY      One of: ['segments', 'bases'], which determines the
                         definition of overlap. @segments: The value associated
                         with two features overlapping will be 1 if they
                         overlap, and 0 otherwise. @bases: The value associated
                         with two features overlapping will be number of base
                         pairs which they overlap. [default: segments]
     --midpoint-only=MIDPOINT
                         For the specified file (1, 2, or both), use onlythe
                         midpoint of each feature instead of the entire width.
     -m MIN_OVERLAP, --min-overlap=MIN_OVERLAP
                         The minimum number of base pairs that two features
                         must overlap for them to be classified as overlapping.
                         This integer can be either positive (features overlap
                         only if they share at least this many bases) or
                         negative (features overlap if there are no more than
                         this many bases between them). Both a negative min-
                         overlap and --by=bases cannot be specified together.
                         [default: 1]
     --min-overlap-fraction=MIN_OVERLAP_FRACTION
                         The minimum fraction of the base pairs in the subject
                         feature that overlap with the query feature in order
                         to be counted as overlapping. Overrides--min-overlap.

   Files:
     --mnemonic-file=MNEMONICFILENAME
                         If specified, labels will be shown using mnemonics
                         found in this file
     -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: overlap]

   GSC Options:
     --region-file=REGIONFILENAME
                         If specified, this file will be used to calculate
                         overlap significance using GSC. This must be a BED
                         file
     -s SAMPLES, --samples=SAMPLES
                         The number of samples for GSC to use to estimate the
                         significance of the overlap [default: 1000]
     --region-fraction=REGION_FRACTION
                         The region_fraction tu use with GSC [default: 0.5]
     --subregion-fraction=SUBREGION_FRACTION
                         The subregion_fraction tu use with GSC [default: 0.5]

:mod:`.signal_distribution`
---------------------------

.. module:: segtools.signal_distribution

Provides command-line and package entry points for analyzing the signal
distribution over tracks and labels.

.. program:: segtools-signal-distribution

Command-line usage summary
..........................

::

 Usage: segtools-signal-distribution [OPTIONS] BEDFILE GENOMEDATADIR

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit

   Flags:
     --clobber           Overwrite existing output files if the specified
                         directory already exists.
     --quick             Compute values only for one chromosome.
     --replot            Load data from output tab files and regenerate plots
                         instead of recomputing data
     --noplot            Do not generate plots
     --group-labels      Group track distributions over all labels. BEDFILE
                         will be ignored
     --ecdf              Plot empiracle cumulative density inside each panel
                         instead of a normal histogram (turns off log-y)
     --calc-ranges       Calculate ranges for distribution plots from
                         segmentation data (slower) instead of using whole
                         genome data (default).

   Histogram options:
     -n NUM_BINS, --num-bins=NUM_BINS
                         Number of bins for signal distribution [default: 100]
     --min-value=MIN_VALUE
                         Minimum signal track value used in binning (overrides
                         min from --calc-ranges) (values below will be ignored)
     --max-value=MAX_VALUE
                         Maximum signal track value used in binning (overrides
                         max from --calc-ranges) (values above will be ignored)

   I/O options:
     --mnemonic-file=MNEMONICFILENAME
                         If specified, labels will be shown using mnemonics
                         found in this file
     -o OUTDIR, --outdir=OUTDIR
                         File output directory (will be created if it does not
                         exist) [default: signal_distribution]  
