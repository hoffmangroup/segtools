======================
Segtools documentation
======================
:Author: Orion Buske <orion.buske at gmail dot com>
:Organization: University of Washington
:Address: Department of Genome Sciences, PO Box 355065, Seattle, WA 98195-5065, United States of America
:Copyright: 2009 Orion Buske

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
   
   > wget http://http://encodestatistics.org/svn/segmentation_validation/segtools/trunk/install.py
   > python install.py

.. _script: http://http://encodestatistics.org/svn/segmentation_validation/segtools/trunk/install.py


Usage
=====
All segtools modules requires, at the very least, a segmentation in 
``BED`` format_. Some modules also require additional files, such as 
genomic feature files to compare with the segmentation. Modules are most
easily used from the command line interface, but all can be loaded and 
run straight from python if desired.

.. _format: http://genome.ucsc.edu/FAQ/FAQformat#format1

The following modules are currently available::
    
    - feature-aggregation_: Analyzes the relative occurrance of each segment
    label around the provided genomic features.
    - label-transition: Analyses the transitions between segment labels and
    the structure of their interaction.
    - length-distribution_: Analyzes the distribution of segment lengths
    and their coverage of the genome for each segment label.
    - signal-distribution_: Analyzes the distribution of genomic signal 
    tracks for each segment label.
    - nucleotide-frequency_: Analyzes the frequencies of nucleotides and
    dinucleotides in segments of each label.
    - overlap_: Analyzes the frequency with which each segment label overlaps
    features of various types.
    - *html-report_:* Combines the output of the other modules and generates
    an html report for viewing.

These modules can be run from the command line by prepending ``segtools-``
to the module name, such as::

   > segtools-length-distribution --help

Each module generates tab-delimited (``tab``) data files, image files 
(in ``png`` and ``pdf`` format and in normal, thumbnail, and 
slide layouts), and partial HTML (``div``) files.


Technical description
---------------------

blah?


Modules
=======


feature-aggregation
-------------------

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


html-report
-------------------

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


label-transition
----------------

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

length-distribution
-------------------

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




nucleotide-frequency
-------------------

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


overlap
-------------------

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

signal-distribution
-------------------

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


Segway generates a model (``segway.str``) and initial parameters
(``input.master``) appropriate to a dataset using the GMTKL
specification language and the GMTK master parameter file format. Both
of these are described more fully in the GMTK documentation (cite),
and the default structure and starting parameters are described more
fully in the Segway article.

You can tell Segway just to generate these files and not to perform
any inference using the ``--dry-run`` option.

Using :option:`--num-starts`\=\ *starts* will generate multiple copies of the
``input.master`` file, named ``input.0.master``, ``input.1.master``,
and so on, with different randomly picked initial parameters. You may
substitute your own ``input.master`` files but I recommend starting
with a Segway-generated template. This will help avoid some common
pitfalls. In particular, if you are going to perform training on your
model, you must ensure that the ``input.master`` file retains the same
``#ifdef`` structure for parameters you wish to train. Otherwise, the
values discovered after one round of training will not be used in
subsequent rounds, or in the identify or posterior stages.

You can use the :option:`--num-labels`\=\ *labels* option to specify the
number of segment labels to use in the model (default 2). You can set
this to a single number or a range with Python slice notation. For
example, ``--num-labels=5:20:5`` will result in 5, 10, and 15 labels
being tried. If you specify :option:`--num-starts`\=\ *starts*, then
there will be *starts* different threads for each of the *labels*
labels tried.

Segway allows multiple models of the values of an observation track
using three different probability distributions: the normal
distribution, the gamma distribution, and a multinomial distribution,
where the nominal output classes each map to a different bin of the
numerical data.

XXX cleanup duplication

The model may be generated using a normal distribution for continuous
observed tracks (``--distribution=norm``, the default), a normal
distribution on asinh-transformed data
(``--distribution=asinh_norm``), or a gamma distribution
(``--distribution=gamma``). The ideal methodology for setting gamma
parameter values is less well-understood, and it also requires an
unreleased version of GMTK. I recommend the use of ``asinh_norm`` in
most cases.

For gamma distributions, Segway generates initial parameters by
converting mean~$\mu$ and variance~$\sigma^2$ to shape~$k$ and
scale~$\theta$ using the equations~$\mu = k \theta$ and~$\sigma^2 = k
\theta^2$.

XXX add arcsinh_normal similar to log

You may specify a subset of tracks using the ``--trackname`` option
which may be repeated. For example::

    segway --trackname dnasei --trackname h3k36me3

will include the two tracks ``dnasei`` and ``h3k36me3`` and no others.

It is very important that you always specify the same ``--trackname``
options at all stages in the Segway workflow. There is also a special
track name, ``dinucleotide``. When you specify
``--trackname=dinucleotide``, Segway will create a track containing
the dinucleotide that starts at a particular position. This can help
in modeling CpG or G+C bias.

Segment length constraints
==========================

The XXX option allows specification of minimum and maximum segment
lengths for various labels. XXX include sample of table

also a way to add a soft prior on XXX cover --prior-strength
XXX default is XXXcomp, this can't be changed at the moment. E-mail
Michael if you need it to be changable.

Distributed computing
=====================
Segway can currently perform training and identification tasks only
using a cluster controllable with the DRMAA (cite) interface. I have
only tested it against Sun Grid Engine, but it should be possible to
work with other DRMAA-compatible distriuted computing systems, such as
Platform LSF, PBS, Condor, (XXXcomp add others). If you are interested
in using one of these systems, please contact me so we can correct all
the fine details. A standalone version is planned.

Training
========
Most users will generate the model at training time, but to specify
your own model there are the ``--structure=<filename>`` and
``--input-master=<filename>`` options.

Training can be a time-consuming process. You may wish to train only
on a subset of your data. To facilitate this, there is an
``--include-regions=<filename>`` option which specifies a BED file
containing a list of regions to limit to. For example, the ENCODE Data
Coordination Center at University of Califronia Santa Cruz keeps the
coordinates of the ENCODE pilot regions in this format at XXXcomp. For
human whole-genome studies, these regions have nice properties since
they mark 1 percent of the genome, and were carefully picked to
include a variety of different gene densities, and a number of more
limited studies provide data just for these regions. There is a file
containing only nine of these regions at XXXcomp(make it), which
covers 0.15% of the human genome, and is useful for training.

Memory usage
============

XXX describe new regime

XXX other sections of workflow

XXX other sections of technical description

XXX add section on all other options

Command-line usage summary
==========================

XXX cover all of these options.

::

  Usage: segway [OPTION]... GENOMEDATADIR
  
  Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit
  
    Data subset:
      -t TRACK, --track=TRACK
                          append TRACK to list of tracks to use (default all)
      --include-coords=FILE
                          limit to genomic coordinates in FILE
      --exclude-coords=FILE
                          filter out genomic coordinates in FILE
  
    Model files:
      -i FILE, --input-master=FILE
                          use or create input master in FILE
      -s FILE, --structure=FILE
                          use or create structure in FILE
      -p FILE, --trainable-params=FILE
                          use or create trainable parameters in FILE
      --dont-train=FILE   use FILE as list of parameters not to train
      --seg-table=FILE    load segment hyperparameters from FILE
      --semisupervised=FILE
                          semisupervised segmentation with labels in FILE
  
    Output files:
      -b FILE, --bed=FILE
                          create bed track in FILE
  
    Intermediate files:
      -o DIR, --observations=DIR
                          use or create observations in DIR
      -d DIR, --directory=DIR
                          create all other files in DIR
  
    Variables:
      -D DIST, --distribution=DIST
                          use DIST distribution
      -r NUM, --random-starts=NUM
                          randomize start parameters NUM times (default 1)
      -N SLICE, --num-segs=SLICE
                          make SLICE segment classes (default 2)
      --prior-strength=RATIO
                          use RATIO times the number of data counts as the
                          number of pseudocounts for the segment length prior
                          (default 0)
      -m PROGRESSION, --mem-usage=PROGRESSION
                          try each float in PROGRESSION as the number of
                          gibibytes of memory to allocate in turn (default
                          2,3,4,6,8,10,12,14,15)
      -v NUM, --verbosity=NUM
                          show messages with verbosity NUM
      --drm-opt=OPT       specify an option to be passed to the distributed
                          resource manager
  
    Flags:
      -c, --clobber       delete any preexisting files
      -T, --no-train      do not train model
      -I, --no-identify   do not identify segments
      -P, --no-posterior  do not identify probability of segments
      -k, --keep-going    keep going in some threads even when you have errors
                          in another
      -n, --dry-run       write all files, but do not run any executables
      -S, --split-sequences
                          split up sequences that are too large to fit into
                          memory
  
  

Python interface
================
I have designed Segway such that eventually one may call different
components directly from within Python. To do so, import the following
module:

XXXcomp table here (from the setup.py)

You can then call the appropriate module through its ``main()``
function with the same arguments you would use at the command line.
For example::

  from segway import run

  GENOMEDATA_DIRNAME = "genomedata"

  run.main("--no-identify", GENOMEDATA_DIRNAME)

XXX describe runner.fromoptions() interface

All other interfaces (the ones that do not use a ``main()`` function)
to Segway code are undocumented and should not be used. If you do use
them, know that the API may change at any time without notice.
