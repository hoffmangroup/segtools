#!/usr/bin/env python
from __future__ import division, with_statement

"""
This module aggregates segmentation data around features, generating
a histogram for each segmentation label that shows the frequency of
observing that label at that position relative the the feature.

If using gene mode, the input file should have features with names:
exon, start_codon, CDS
as provided by exporting UCSC gene data in GTF format.
"""

# A package-unique, descriptive module name used for filenames, etc
MODULE="feature_aggregation"

import os
import sys

from collections import defaultdict
from functools import partial
from numpy import arange, array, empty, linspace, round

from . import (Annotation, log, Segmentation, DTYPE_SEGMENT_START,
               DTYPE_SEGMENT_END, DTYPE_SEGMENT_KEY, DTYPE_STRAND,
               die, RInterface, open_transcript, add_common_options)
from .common import (fill_array, get_ordered_labels, make_tabfilename,
                     map_segment_label, setup_directory, tab_saver)
from .html import save_html_div
from .version import __version__

NAMEBASE = "%s" % MODULE
NAMEBASE_SPLICING=os.extsep.join([NAMEBASE, "splicing"])
NAMEBASE_TRANSLATION=os.extsep.join([NAMEBASE, "translation"])

STATIC_FIELDNAMES = ["group", "component", "offset"]

HTML_TITLE_BASE = "Feature aggregation"

# REGION_MODE, GENE_MODE mean do quantile bins
POINT_MODE = "point"
REGION_MODE = "region"
GENE_MODE = "gene"
MODES = [POINT_MODE, REGION_MODE, GENE_MODE]
DEFAULT_MODE = POINT_MODE

# Component names
#   A single %f/%d/%s will automatically be substituted with the average number
#   of bases in that component over all features
FLANK_5P = "5' flanking: %d bp"
FLANK_3P = "3' flanking: %d bp"
FLANK_COMPONENTS = [FLANK_5P, FLANK_3P]
POINT_COMPONENTS = FLANK_COMPONENTS

REGION_COMPONENT = "internal (%d bp)"
REGION_COMPONENTS = [FLANK_5P, REGION_COMPONENT, FLANK_3P]

INITIAL_EXON = "initial exon (%d bp)"
TERMINAL_EXON = "terminal exon (%d bp)"
SPLICE_COMPONENTS = [
    INITIAL_EXON,
    "initial intron (%d bp)",
    "internal exons (%d bp)",
    "internal introns (%d bp)",
    TERMINAL_EXON,
    "terminal intron (%d bp)"]
CODING_COMPONENTS = [
    "initial 5' UTR (%d bp)",
    "5' UTR introns (%d bp)",
    "internal 5' UTR (%d bp)",
    "terminal 5' UTR (%d bp)",
    "initial CDS (%d bp)",
    "terminal CDS (%d bp)",
    "initial 3' UTR (%d bp)",
    "internal 3' UTR (%d bp)",
    "3' UTR introns (%d bp)",
    "terminal 3' UTR (%d bp)"]

EXON_COMPONENTS = [FLANK_5P] + SPLICE_COMPONENTS + [FLANK_3P]
GENE_COMPONENTS = EXON_COMPONENTS + CODING_COMPONENTS

# For parsing GTF file
EXON_FEATURE = "exon"
CDS_FEATURE = "CDS"
MIN_EXONS = 1
MIN_CDSS = 0

# Default parameter values
FLANK_BASES = 500
REGION_BINS = 50
INTRON_BINS = 50
EXON_BINS = 25

R = RInterface(["common.R", "aggregation.R"])

class GeneAnnotation(Annotation):
    """Annotation class for gene-type annotation"""

    def _from_file(self, filename, verbose=True):
        """Load the gtf file in terms of idealized gene features.

        Parses the given feature file and replaces feature names from:
          CDS
          exon
          (other features ignored)
        with feature components from:
          GENE_COMPONENTS

        Field 8 of the feature file must be a list of semi-colon delimited
        properties and gene_id must be the first one of them.

        All features must be stranded.

        """
        from .common import inverse_dict

        log("  Parsing GTF file...", verbose, end="")
        format = self._get_file_format(filename)
        if format != "gtf":
            raise self.FormatError("ANNOTATION must be in GTF format to"
                                   " aggregate with --mode=gene")

        # Start with establishing a dict:
        #   gene_id -> dict(transcript_id -> (chrom, start, end, name, strand))
        #   start is 0-indexed, end is non-inclusive (BED)
        gene_dict = defaultdict(partial(defaultdict, list))
        for row in self._iter_rows(filename, verbose=verbose):
            if row['strand'] not in {"+", "-"}:
                raise self.FormatError("ANNOTATION must specify a strand"
                                       " for every entry to aggregate with"
                                       " --mode=gene")

            # Skip GTF lines where the feature field == gene
            # Prevent gene lines to be selected when
            # calculating the longest transcript
            if row['name'] != "gene":
                # Add feature to dict
                gene_dict[row['gene_id']][row['transcript_id']].append(row)

        log(" done", verbose)
        log("  Mapping onto gene model...", verbose, end="")

        # Eventually create feature dict: chrom -> list(features)
        data = defaultdict(list)
        label_dict = {}
        n_genes = 0
        n_label_segments = {}
        n_label_bases = {}
        for gene_id, transcript_dict in gene_dict.iteritems():
            # Select only longest transcript
            transcript_id, transcript = \
                self._get_longest_transcript(transcript_dict)
            gene_model = self._transcript_to_gene_features(transcript)

            if gene_model is not None:
                n_genes += 1
                # Add features to a normal, chrom-based feature dict
                for gene_part in gene_model:
                    # label is the component of that gene part
                    chrom, start, end, strand, label = gene_part
                    try:
                        label_key = label_dict[label]
                    except KeyError:
                        label_key = len(label_dict)
                        label_dict[label] = label_key
                        n_label_segments[label] = 0
                        n_label_bases[label] = 0

                    feature = (start, end, label_key, strand)
                    data[chrom].append(feature)
                    n_label_segments[label] += 1
                    n_label_bases[label] += end - start
#            raw_input()

        labels = inverse_dict(label_dict)
        dtype = [('start', DTYPE_SEGMENT_START),
                 ('end', DTYPE_SEGMENT_END),
                 ('key', DTYPE_SEGMENT_KEY),
                 ('strand', DTYPE_STRAND)]

        chromosomes = dict((chrom, array(features, dtype=dtype))
                           for chrom, features in data.iteritems())

        log(" done", verbose)
        log("  Sorting...", verbose, end="")

        # Sort features by ascending start
        for chrom_features in chromosomes.itervalues():
            chrom_features.sort(order=['start'])

        log(" done", verbose)

        self.filename = filename
        self.chromosomes = chromosomes
        self._labels = labels
        self._n_label_segments = n_label_segments
        self._n_label_bases = n_label_bases
        self._n_genes = n_genes
        self._inv_labels = label_dict

    def num_genes(self):
        return self._n_genes

    def _preprocess_entries(self, rows):
        """Convert list of gene feature entries to processed gene data.

        If multiple transcripts are found, only the longest is taken.

        :param entries: each entry is a dict
                        and all entries should correspond to the same gene_id
        :returns: chrom, strand, list(exon), list(cds)

        """
        # Remove list of gene features and preprocess
        #   (extracting constant strand, chrom)
        # Returns gene info and list of exons and CDSs
        gene_strand = None
        gene_chrom = None
        exons = []
        cdss = []

        for row in rows:
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            strand = row['strand']
            name = row['name']

            # Ensure strand and chrom match rest for this gene
            if gene_strand is None:
                gene_strand = strand
            elif gene_strand != strand:
                die("Found gene features on more than one strand: [%s, %s]\n%s"
                    % (gene_strand, strand, row))

            if gene_chrom is None:
                gene_chrom = chrom
            elif gene_chrom != chrom:
                die("Found gene features on more than one chromosome:"
                    " [%s, %s]\n%s" % (gene_chrom, chrom, row))

            partial_entry = (start, end)
            if name == EXON_FEATURE:
                exons.append(partial_entry)
            elif name == CDS_FEATURE:
                cdss.append(partial_entry)
        return gene_chrom, gene_strand, exons, cdss

    def _transcript_to_gene_features(self, rows):
        """Returns a set of gene-model entries based upon the given entries.

        Interprets the given set of entries in terms of an idealized gene model
        with 5' and 3' UTRs, initial, internal, and terminal exons and introns.

        :param entries: see preprocess_entries for a description
        :returns: a new set of entries based upon this gene model
                  or None if the entries don't fit the model.
        """

        chrom, strand, exons, cdss = self._preprocess_entries(rows)

        # Ignore genes without any exons
        if not exons:
            return None

        exons.sort()
        cdss.sort()

        # Create introns between every pair of now-sorted exons
        introns = []
        prev_end = None
        for exon_start, exon_end in exons:
            if prev_end is not None:
                introns.append((prev_end, exon_start))
            prev_end = exon_end

        # Flip directions for '-' strand
        if strand == "-":
            exons.reverse()
            introns.reverse()
            cdss.reverse()

        features = []
        # Add gene features to new list, renaming components
        #   according to an idealized gene model
        # Assumes exons and cdss in sorted order (ensured by sort() above)

        # Binds local variables: chrom, strand
        def add_feature(component, start, end):
            """Makes tuple feature for feature list"""
            features.append((chrom, start, end, strand, component))  # tuple

        # first exon
        add_feature(SPLICE_COMPONENTS[0], *exons[0])
        # first intron
        if len(introns) > 0:
            add_feature(SPLICE_COMPONENTS[1], *introns[0])
        # internal exons
        for exon in exons[1:-1]:
            add_feature(SPLICE_COMPONENTS[2], *exon)
        # internal introns
        for intron in introns[1:-1]:
            add_feature(SPLICE_COMPONENTS[3], *intron)
        # last intron
        if len(introns) > 0:
            add_feature(SPLICE_COMPONENTS[5], *introns[-1])
        # last exon
        add_feature(SPLICE_COMPONENTS[4], *exons[-1])


        # If there were no CDS's, then there is no point in dealing with UTR/CDS
        #   regions, so skip all that.
        if len(cdss) == 0:
            return features

        # Binds local variable strand to strand-correct
        def upstream(x, y):
            """Returns true if x starts before (is 5') of y, false otherwise"""
            return (strand == "+" and x[0] < y[0]) or \
                (strand == "-" and x[1] > y[1])

        # Go through exons and introns and pick out those that are in UTR
        # Exons are a little trickier because the UTR region has to be trimmed
        #   to the CDS if they overlap
        UTR5p_exons = []
        UTR5p_introns = []
        UTR3p_exons = []
        UTR3p_introns = []
        first_cds = cdss[0]
        last_cds = cdss[-1]
        first_cds_start, first_cds_end = first_cds
        last_cds_start, last_cds_end, = last_cds
        for exon in exons:
            exon_start, exon_end = exon
            # See if exon is in UTR relative to CDSs, and trim to UTR fragment
            utr_list = None
            if strand == "+":
                if exon_start < first_cds_start:
                    exon_end = min(exon_end, first_cds_start)
                    utr_list = UTR5p_exons
                elif exon_end > last_cds_end:
                    exon_start = max(exon_start, last_cds_end)
                    utr_list = UTR3p_exons
            else:  # strand == "-":
                if exon_end > first_cds_end:
                    exon_start = max(exon_start, first_cds_end)
                    utr_list = UTR5p_exons
                elif exon_start < last_cds_start:
                    exon_end = min(exon_end, last_cds_start)
                    utr_list = UTR3p_exons
            # Add if exon was in UTR and has positive length
            if utr_list is not None and exon_start < exon_end:
                utr_list.append((exon_start, exon_end))  # Append tuple

        for intron in introns:
            if upstream(intron, first_cds):
                UTR5p_introns.append(intron)
            elif upstream(last_cds, intron):
                UTR3p_introns.append(intron)


        # Go through the list of UTR elements and add them to the feature list
        #   with appropriate component names
        # initial 5' UTR
        if len(UTR5p_exons) > 0:
            add_feature(CODING_COMPONENTS[0], *UTR5p_exons[0])
        # 5' UTR introns
        for intron in UTR5p_introns:
            add_feature(CODING_COMPONENTS[1], *intron)
        # internal 5' UTR
        for exon in UTR5p_exons[1:-1]:
            add_feature(CODING_COMPONENTS[2], *exon)
        # terminal 5' UTR
        if len(UTR5p_exons) > 0:
            add_feature(CODING_COMPONENTS[3], *UTR5p_exons[-1])
        # first CDS
        add_feature(CODING_COMPONENTS[4], *first_cds)
        # last CDS
        add_feature(CODING_COMPONENTS[5], *last_cds)
        # initial 3' UTR
        if len(UTR3p_exons) > 0:
            add_feature(CODING_COMPONENTS[6], *UTR3p_exons[0])
        # internal 3' UTR
        for exon in UTR3p_exons[1:-1]:
            add_feature(CODING_COMPONENTS[7], *exon)
        # 3' UTR introns
        for intron in UTR3p_introns:
            add_feature(CODING_COMPONENTS[8], *intron)
        # terminal 3' UTR
        if len(UTR3p_exons) > 0:
            add_feature(CODING_COMPONENTS[9], *UTR3p_exons[-1])

#         for row in rows:
#             log("%s\t%d\t%d\t%s" % (row['chrom'], row['start'],
#                                     row['end'], row['name']))
#         for feature in features:
#             log("\t%s" % str(feature))
#         raw_input()
        return features

    @staticmethod
    def _get_longest_transcript(transcript_dict):
        """Return the longest transcript in the dict.

        :type transcript_dict: transcript_id -> entries
        :returns: tuple of (transcript_id, entries) for the longest transcript
                  or None if there were no transcripts.
        """

        def get_transcript_length(rows):
            start = min(row['start'] for row in rows)
            end = max(row['end'] for row in rows)

            return end - start

        # XXX: this could probably be simplified vastly using max()
        # and a generator expression; the semantics would change
        # slightly in case of ties; you would get the transcript with
        # the highest ID instead of the first one of that length.
        # check to see if this is OK
        max_length = 0
        res = None
        for id, rows in transcript_dict.iteritems():
            length = get_transcript_length(rows)
            if length > max_length:
                max_length = length
                res = (id, rows)
        return res


def print_bed_from_gene_component(features, component="terminal exon"):
    with open("%s.bed" % component, "w") as ofp:
        for chrom, chrom_features in features.iteritems():
            for feature in chrom_features:
                if feature["name"].startswith(component):
                    print >>ofp, "%s\t%s\t%s\t%s" % (chrom,
                                                     feature["start"],
                                                     feature["end"],
                                                     feature["name"])


def calc_feature_windows(feature, labels, mode, component_bins):
    """Calculate where the bin windows for this feature begin; Spread
    internal bins throughout feature

    Return a list of tuples: (component, window)

    component: str (usually from *_COMPONENT/*_COMPONENTS constants)
    window: numpy.ndarray(dtype=int) of base indices
    """
    name = labels[feature['key']]
    start = feature['start']
    end = feature['end']
    try:
        strand = feature['strand']
    except (IndexError, KeyError):
        assert mode != GENE_MODE
        strand = "."

    length = end - start
    assert length >= 0

    # Include flanks by default
    num_5p_bins = component_bins[FLANK_5P]
    num_3p_bins = component_bins[FLANK_3P]

    # remove flanks if gene mode and not on one end or other
    if mode == GENE_MODE:
        assert name in component_bins
        component = name
        if component != INITIAL_EXON:
            # Not initial exon, so no 5' flank
            num_5p_bins = 0

        if component != TERMINAL_EXON:
            # Not terminal exon, so no 3' flank
            num_3p_bins = 0
    elif mode == REGION_MODE:
        component = REGION_COMPONENT
    else:
        assert mode == POINT_MODE
        component = None
        num_internal_bins = 0

    if component is not None:
        # default value of 0
        num_internal_bins = component_bins.get(component, 0)

    if num_internal_bins > length:
        # don't calculate internal bins in this case
        #print >> sys.stderr, "Warning: %d %s bins > %d bases" % \
        #    (num_internal_bins, component, length)
        num_internal_bins = 0

    assert num_internal_bins >= 0

    # calculate internal bin locations
    unrounded_bins = linspace(start, end, num_internal_bins, endpoint=False)
    internal_bins = round(unrounded_bins).astype(int)

    # Calculate flanking base locations
    # XXX: what happens when the flanking regions are off the chromosome?
    if strand == "-":
        bins_5p = arange(end, end + num_5p_bins, dtype="int")[::-1]
        internal_bins = internal_bins[::-1]  # reverse
        bins_3p = arange(start - num_3p_bins, start, dtype="int")[::-1]
    else:
        bins_5p = arange(start - num_5p_bins, start, dtype="int")
        bins_3p = arange(end, end + num_3p_bins, dtype="int")

    res = []
    if len(bins_5p) > 0:
        res.append((FLANK_5P, bins_5p))
    if len(internal_bins) > 0:
        res.append((component, internal_bins))
    if len(bins_3p) > 0:
        res.append((FLANK_3P, bins_3p))

    #for key, val in locals().iteritems():
    #    print >>sys.stderr, "%s : %s" % (key, val)
    #sys.exit(1)

    return res

## Accepts data from dict: chr -> dict {"start", "end", "strand"})
##   zero-based start and end (exclusive) indices
## If components is:
##   [] or None, aggregation will be just over the flanking regions
##   length 1, aggregation will be in "region" mode
##   length > 1, aggregation will be over each component separately, with
##     flanking regions before the 1st and after the last component in the list
##     Each feature's component entry must match one of these exactly
## component_bins is a dict: component -> number of bins for that component
## If not by_groups: all groups found are treated as only group in groups
## Returns:
##   groups: a list of the groups aggregated over
##   counts: dict(label_key -> dict(feature -> dict(component -> histogram)))
def calc_aggregation(segmentation, features, mode, groups, components,
                     component_bins, quick=False, by_groups=False,
                     verbose=True):
    if by_groups:
        assert len(groups) > 0
    else:
        assert len(groups) == 1

    labels = segmentation.labels

    if verbose:
        log("\tGroups: %s" % groups)
        log("\tComponents and bins:")
        for component in components:
            log("\t\t", end="")
            log("%s; %d bins" % (component, component_bins[component]))

    # counts is an accumulator
    # dict:
    #   key: feature_group
    #   value: dict:
    #      key: component_name
    #      value: numpy.ndarray histogram [bin, label_key]
    counts = dict([(group,
                    dict([(component, fill_array(0, (bins, len(labels))))
                          for component, bins in component_bins.iteritems()]))
                   for group in groups])

    counted_features = 0

    for chrom, segments in segmentation.chromosomes.iteritems():
        log("\t%s" % chrom, verbose)

        try:
            chrom_features = features.chromosomes[chrom]
        except KeyError:
            continue

        if len(chrom_features) == 0:
            continue

        # Get bounds on segmentation for chr (since segments are sorted)
        segmentation_start = segments['start'][0]
        segmentation_end = segments['end'][-1]

        # Map entire chromosome's segments into an array of label_keys
        # XXXopt: this line takes up to 3.2 GB of memory
        segment_map, sentinel = \
            map_segment_label(segments, (segmentation_start, segmentation_end))

        # For each feature, tally segments in window
        for feature in chrom_features:
            # XXXopt: Process every feature, since we don't know how wide the
            # aggregation window is.
            feature_counted = False
            if by_groups and mode != GENE_MODE:
                group = features.labels[feature['key']]
            else:
                group = groups[0]

            group_counts = counts[group]
            # dict:
            #      key: component_name
            #      value: numpy.ndarray histogram [bin, label_key]

            # Spread internal bins throughout feature if not POINT_MODE
            component_windows = calc_feature_windows(feature, features.labels,
                                                     mode, component_bins)

            # Scan window, tallying segments observed
            for component, window in component_windows:
                component_counts = group_counts[component]
                #log("component: %r" % component)
                #log("component_counts.shape: %r" % str(component_counts.shape))

                map_indices = window - segmentation_start
                count_indices = arange(0, len(window))

                # Lookup all label_keys within segmentation at once
                keep = ((window >= segmentation_start) &
                        (window < segmentation_end))
                label_keys = segment_map[map_indices[keep]]
                if len(label_keys) == 0:
                    continue
                count_indices = count_indices[keep]

                # Find indices where label_keys are set
                keep = (label_keys != sentinel).nonzero()[0]
                label_keys = label_keys[keep]
                if len(label_keys) == 0:
                    continue
                count_indices = count_indices[keep]

                #raw_input()
                component_counts[count_indices, label_keys] += 1
                if len(keep) > 0 and not feature_counted:
                    counted_features += 1
                    feature_counted = True

        if quick:
            break

    return (counts, counted_features)

def make_row(labels, row_data):
    values = {}
    for label_key, label in labels.iteritems():
        values[label] = row_data[label_key]

    return values

## Saves the data to a tab file
def save_tab(segmentation, features, counts, components, component_bins,
             counted_features, dirpath, mode, verbose=True,
             clobber=False, namebase=NAMEBASE):
    metadata = {"num_features": counted_features}
    # Add metadata to tell R how to display gene model.
    if mode == GENE_MODE:
        metadata["spacers"] = len(EXON_COMPONENTS)

    label_keys, labels = get_ordered_labels(segmentation.labels)
    for label_key in label_keys:
        label = labels[label_key]
        assert label not in STATIC_FIELDNAMES
        assert label not in metadata
        metadata[label] = segmentation.num_label_bases(label)

    log("Saving metadata: %r" % metadata, verbose)

    fieldnames = STATIC_FIELDNAMES + [labels[label_key]
                                      for label_key in label_keys]
    with tab_saver(dirpath, namebase, fieldnames=fieldnames, metadata=metadata,
                   clobber=clobber, verbose=verbose) as saver:
        # Compute the average length for features in each component,
        #   merging across groups.
        avg_component_lengths = {}
        for component in components:
            avg_len = None
            if component in FLANK_COMPONENTS:
                avg_len = component_bins[component]
            else:
                sum = 0
                count = 0
                for chrom_features in features.chromosomes.itervalues():
                    if component == REGION_COMPONENT:
                        # Merge groups
                        sum += (chrom_features['end'] -
                                chrom_features['start']).sum()
                        count += len(chrom_features)
                    elif component in GENE_COMPONENTS:
                        # Separate by component (key column)
                        try:
                            label_key = features.label_key(component)
                        except KeyError:
                            continue

                        keep = chrom_features['key'] == label_key
                        sum += (chrom_features['end'][keep] -
                                chrom_features['start'][keep]).sum()
                        count += keep.sum()

                if sum > 0:
                    avg_len = sum / count

            avg_component_lengths[component] = avg_len

        for group in counts:
            for component in components:
                hist = counts[group][component]
                # Try to substitute average lengths into component names
                component_name = component
                avg_len = avg_component_lengths[component]
                if avg_len is not None:
                    try:
                        component_name = component % avg_len
                    except TypeError:
                        pass

                if component == FLANK_5P:
                    offsets = range(-len(hist), 0)
                else:
                    offsets = range(0, len(hist))

                for offset, row_data in zip(offsets, hist):
                    row = make_row(labels, row_data)
                    row["group"] = group
                    row["component"] = component_name
                    row["offset"] = offset
                    saver.writerow(row)

## Plots aggregation data from tab file
def save_plot(dirpath, mode, namebase=NAMEBASE, clobber=False, verbose=True,
              mnemonic_file=None, normalize=False, significance=False,
              transcriptfile=None, ropts=None):
    R.start(verbose=verbose, transcriptfile=transcriptfile)

    tabfilename = make_tabfilename(dirpath, namebase)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    if not mnemonic_file:
        mnemonic_file = ""

    # Plot data in tab file
    if mode == GENE_MODE:
        R.plot("save.gene.aggregations", dirpath, NAMEBASE_SPLICING,
               NAMEBASE_TRANSLATION, tabfilename, mnemonic_file=mnemonic_file,
               normalize=normalize, clobber=clobber, significance=significance,
               ropts=ropts)
    else:
        R.plot("save.aggregation", dirpath, namebase, tabfilename,
               mnemonic_file=mnemonic_file, normalize=normalize,
               clobber=clobber, significance=significance, ropts=ropts)

def save_html(dirpath, featurefilename, mode, clobber=False,
              verbose=True, normalize=False):
    featurebasename = os.path.basename(featurefilename)
    title = "%s (%s)" % (HTML_TITLE_BASE, featurebasename)
    if normalize:
        yaxis = "enrichment"
    else:
        yaxis = "count"

    if mode == GENE_MODE:
        templatefile = "aggregation_gene_div.tmpl"
        extra_namebases = {"splice": NAMEBASE_SPLICING,
                           "trans": NAMEBASE_TRANSLATION}
    else:
        templatefile = "aggregation_div.tmpl"
        extra_namebases = {}

    save_html_div(templatefile, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, featurefilename=featurebasename, mode=mode,
                  title=title, yaxis=yaxis, verbose=verbose,
                  extra_namebases=extra_namebases)

def print_array(arr, tag="", type="%d"):
    if len(arr) > 2:
        fstring = "%%s: \t%s, %s, ..., %s, ..., %s, %s" % tuple([type]*5)
        log(fstring % (tag,
                       arr[0],
                       arr[1],
                       arr[int(arr.shape[0] / 2)],
                       arr[-2],
                       arr[-1]))
    else:
        log(str(arr))

def get_components_and_bins(mode, flank_bases=FLANK_BASES,
                            region_bins=REGION_BINS,
                            intron_bins=INTRON_BINS,
                            exon_bins=EXON_BINS):
    """Get components and number of data points per component"""
    if mode == GENE_MODE:
        components = GENE_COMPONENTS
    elif mode == REGION_MODE:
        components = REGION_COMPONENTS
    else:
        assert mode == POINT_MODE
        components = POINT_COMPONENTS

    component_bins = {}
    for component in components:
        if component in FLANK_COMPONENTS:
            bins = flank_bases
        elif component in GENE_COMPONENTS:
            if "intron" in component:
                bins = intron_bins
            elif "exon" in component or \
                    "UTR" in component or \
                    "CDS" in component:
                # All UTR components that don't contain "intron" are exons
                bins = exon_bins
            else:
                bins = region_bins
        else:
            bins = region_bins

        assert component not in component_bins  # Keep unique
        component_bins[component] = bins

    return components, component_bins

## Package entry point
def validate(bedfilename, featurefilename, dirpath,
             flank_bases=FLANK_BASES, region_bins=REGION_BINS,
             intron_bins=INTRON_BINS, exon_bins=EXON_BINS,
             by_groups=False, mode=DEFAULT_MODE, clobber=False,
             quick=False, replot=False, noplot=False, normalize=False,
             significance=False, mnemonic_file=None, verbose=True,
             ropts=None):

    if not replot:
        setup_directory(dirpath)
        segmentation = Segmentation(bedfilename, verbose=verbose)

        if mode == GENE_MODE:
            features = GeneAnnotation(featurefilename, verbose=verbose)
            num_features = features.num_genes()
            groups = ["genes"]
        else:
            features = Annotation(featurefilename, verbose=verbose)
            num_features = features.num_segments()
            if by_groups:
                groups = features.labels.values()
            else:
                groups = ["features"]

        log("Aggregating over %d features" % num_features, verbose)

        components, component_bins = \
            get_components_and_bins(mode,
                                    flank_bases=flank_bases,
                                    region_bins=region_bins,
                                    intron_bins=intron_bins,
                                    exon_bins=exon_bins)

        res = calc_aggregation(segmentation, features, mode=mode,
                               groups=groups, components=components,
                               component_bins=component_bins, quick=quick,
                               by_groups=by_groups, verbose=verbose)
        counts, counted_features = res

        save_tab(segmentation, features, counts, components, component_bins,
                 counted_features, dirpath, mode, clobber=clobber,
                 verbose=verbose)

    if not noplot:
        with open_transcript(dirpath, MODULE) as transcriptfile:
            save_plot(dirpath, mode, clobber=clobber, verbose=verbose,
                      mnemonic_file=mnemonic_file, normalize=normalize,
                      significance=significance, transcriptfile=transcriptfile,
                      ropts=ropts)

    save_html(dirpath, featurefilename, mode=mode, verbose=verbose,
              clobber=clobber, normalize=normalize)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] SEGMENTATION ANNOTATION"
    description = ("Plot the enrichment of the SEGMENTATION labels relative"
                   " to the position of features in ANNOTATION."
                   " Features can be grouped by the"
                   " 'name'/'feature' column by supplying --groups.")
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Output options")
    add_common_options(group, ['mnemonic_file', 'outdir'], MODULE=MODULE)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet', 'quick', 'replot',
                               'noplot'])

    group.add_option("--groups", action="store_true",
                     dest="by_groups", default=False,
                     help="Separate data into different groups based upon"
                     " ANNOTATION's 'name'/'feature' field"
                     " if --mode=region or --mode=point. Does"
                     " nothing if --mode=gene.")
    group.add_option("--normalize", action="store_true",
                     dest="normalize", default=False,
                     help="Plot the relative frequency of SEGMENTATION labels,"
                     " normalized by the number of segments in that group"
                     " instead of the raw counts"
                     " (normalize over SEGMENTATION labels)")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Main aggregation options")
    group.add_option("--mode", choices=MODES,
                     dest="mode", type="choice", default=DEFAULT_MODE,
                     help="one of: %s. [default: %%default]" % \
                     ", ".join(MODES))
    group.add_option("-f", "--flank-bases", metavar="N",
                     dest="flank_bases", type="int", default=FLANK_BASES,
                     help="Aggregate this many base pairs off each"
                     " end of feature or gene [default: %default]")
    group.add_option("-r", "--region-samples", type="int", metavar="N",
                     dest="region_bins", default=REGION_BINS,
                     help="If --mode=region, aggregate over each internal"
                     " feature by taking this many evenly-spaced samples."
                     " Warning: features with a length < N will be excluded"
                     " from the results [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Gene aggregation options")
    group.add_option("-i", "--intron-samples", type="int", metavar="N",
                     dest="intron_bins", default=INTRON_BINS,
                     help="Aggregate over each intron"
                     " by taking this many evenly-spaced samples"
                     " [default: %default]")
    group.add_option("-e", "--exon-samples", type="int", metavar="N",
                     dest="exon_bins", default=EXON_BINS,
                     help="Aggregate over each exon"
                     " by taking this many evenly-spaced samples"
                     " [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "R options")
    add_common_options(group, ['ropts'])
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) != 2:
        parser.error("Inappropriate number of arguments")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]
    featurefilename = args[1]
    kwargs = dict(options.__dict__)
    outdir = kwargs.pop('outdir')

    validate(bedfilename, featurefilename, outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
