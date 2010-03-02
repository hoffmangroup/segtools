#!/usr/bin/env python
from __future__ import division, with_statement

"""
This module aggregates segmentation data around features, generating
a histogram for each segmentation label that shows the frequency of
observing that label at that position relative the the feature.

If using gene mode, the input file should have features with names:
exon, start_codon, CDS
as provided by exporting UCSC gene data in GFF format.
"""

# A package-unique, descriptive module name used for filenames, etc
# Must be the same as the folder containing this script
MODULE="feature_aggregation"

__version__ = "$Revision$"

import math
import os
import sys

from collections import defaultdict
from functools import partial
from operator import itemgetter
from rpy2.robjects import r, numpy2ri

from .common import die, fill_array, get_ordered_labels, image_saver, load_gff_data, load_segmentation, make_tabfilename, map_mnemonics, map_segments, maybe_gzip_open, r_source, setup_directory, tab_saver

from .html import list2html, save_html_div

NAMEBASE = "%s" % MODULE
STATIC_FIELDNAMES = ["group", "component", "offset"]

HTML_TITLE_BASE = "Feature aggregation"
HTML_TEMPLATE_FILENAME = "aggregation_div.tmpl"

POINT_MODE = "point"
REGION_MODE = "region"
GENE_MODE = "gene"
MODES = [POINT_MODE, REGION_MODE, GENE_MODE]
DEFAULT_MODE = POINT_MODE

# Component names
#   A single %d/%s in the name will automatically be substituted with the number
#   of bins used for that component
FLANK_5P = "5' flanking: %d bp"
FLANK_3P = "3' flanking: %d bp"
FLANK_COMPONENTS = [FLANK_5P, FLANK_3P]
POINT_COMPONENTS = list(FLANK_COMPONENTS)
REGION_COMPONENT = "internal"
REGION_COMPONENTS = [FLANK_5P, REGION_COMPONENT, FLANK_3P]
INITIAL_EXON = "initial exon"
TERMINAL_EXON = "terminal exon"
SPLICE_COMPONENTS = [
    INITIAL_EXON,
    "initial intron",
    "internal exons",
    "internal introns",
    TERMINAL_EXON,
    "terminal intron"]
CODING_COMPONENTS = [
    "initial 5' UTR",
    "5' UTR introns",
    "internal 5' UTR",
    "terminal 5' UTR",
    "initial CDS",
    "terminal CDS",
    "initial 3' UTR",
    "internal 3' UTR",
    "3' UTR introns",
    "terminal 3' UTR"]
EXON_COMPONENTS = [FLANK_5P] + SPLICE_COMPONENTS + [FLANK_3P]
GENE_COMPONENTS = EXON_COMPONENTS + CODING_COMPONENTS

# For parsing GTF file
EXON_COMPONENT = "exon"
CDS_COMPONENT = "CDS"
MIN_EXONS = 1
MIN_CDSS = 0

# Default parameter values
FLANK_BINS = 500
REGION_BINS = 50
INTRON_BINS = 50
EXON_BINS = 25
PNG_BLOCK_SIZE = 200


def start_R():
    r_source("common.R")
    r_source("aggregation.R")

def print_bed_from_gene_component(features, component="terminal exon"):
    with open("%s.bed" % component, "w") as ofp:
        for chrom, chrom_features in features.iteritems():
            for feature in chrom_features:
                if feature["name"].startswith(component):
                    print >> ofp, "%s\t%s\t%s\t%s" % (chrom,
                                                      feature["start"],
                                                      feature["end"],
                                                      feature["name"])

# Given a feature, returns the group it is in
def get_feature_group(feature, mode):
    if mode == GENE_MODE:
        return feature["source"]
    else:
        return feature["name"]


## Returns a dict of component -> window (list of base indices)
def calc_feature_windows(feature, mode, components, component_bins):
    name = feature["name"]
    start = feature["start"]
    end = feature["end"]
    strand = feature["strand"]
    length = end - start
    assert length >= 0

    # Include flanks by default
    num_5p_bins = component_bins[FLANK_5P]
    num_3p_bins = component_bins[FLANK_3P]
    if mode == GENE_MODE:  # by gene
        assert name in components
        component = name
        if component != INITIAL_EXON:
            # Not initial exon, so no 5' flank
            num_5p_bins = 0
        if component != TERMINAL_EXON:
            # Not terminal exon, so no 3' flank
            num_3p_bins = 0
    elif mode == REGION_MODE:  # by region
        component = REGION_COMPONENT
    else:  # by point
        component = None

    try:
        num_internal_bins = component_bins[component]
    except KeyError:
        num_internal_bins = 0

    if num_internal_bins > length:
        #print >> sys.stderr, "Warning: %d %s bins > %d bases" % \
        #    (num_internal_bins, component, length)
        num_internal_bins = 0

    # Calculate internal bin locations
    spacing = length / (num_internal_bins + 1)
    internal_bins = []
    for i in range(0, num_internal_bins):
        # Round to nearest bin
        # XXX: This probably isn't quite right in every case, but it's very close
        internal_bins.append(int(start + ((i + 1)*spacing) + 0.5))

    bins_5p = []
    bins_3p = []
    if strand == "-":
        for i in range(0, num_5p_bins):
            bins_5p.append(end + i)
        for i in range(0, num_3p_bins):
            bins_3p.append(start - num_3p_bins + i)
    else:
        for i in range(0, num_5p_bins):
            bins_5p.append(start - num_5p_bins + i)
        for i in range(0, num_3p_bins):
            bins_3p.append(end + i)

    if strand == "-":  # Flip for negative strand
        bins_5p = reversed(bins_5p)
        internal_bins = reversed(internal_bins)
        bins_3p = reversed(bins_3p)

    windows = {}
    if len(bins_5p) > 0:
        windows[FLANK_5P] = bins_5p
    if len(internal_bins) > 0 and component is not None:
        windows[component] = internal_bins
    if len(bins_3p) > 0:
        windows[FLANK_3P] = bins_3p

    #for key, val in locals().iteritems():
    #    print >>sys.stderr, "%s : %s" % (key, val)
    #sys.exit(1)

    return windows

# Returns a dict: component -> number of bins for that component
# and a list of the new components (with flanking regions)
def get_component_bins(components=[], flank_bins=FLANK_BINS,
                       region_bins=REGION_BINS, intron_bins=INTRON_BINS,
                       exon_bins=EXON_BINS):
    component_bins = {}
    for component in components:
        if component in FLANK_COMPONENTS:
            bins = flank_bins
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

        component_bins[component] = bins

    return component_bins


## Given a list of features, returns:
##   A list of the values that occur for that field
def feature_field_values(features, field="name"):
    vals = []
    val_set = set()
    for chrom_features in features.itervalues():
        for feature in chrom_features:
            val = feature[field]
            if val not in val_set:
                val_set.add(val)
                vals.append(val)

    return vals


def preprocess_entries(entries):
    """Convert list of gene feature entries to processed gene data.

    If multiple transcripts are found, only the longest is taken.

    :param entries: each entry is (chrom, start, end, strand, component, source)
                    and all entries should correspond to the same gene_id
    :type entries: list of 6-tuples
    :returns: chrom, strand, source, list(exons), list(cds's)

    """
    # Remove list of gene features and preprocess
    #   (extracting constant strand, chrom, source)
    # Returns gene info and list of exons and CDSs
    gene_strand = None
    gene_chrom = None
    gene_source = None
    exons = []
    cdss = []

    for entry in entries:
        #print entry[4], "\t", entry[1:3], entry[3]
        chrom, start, end, strand, component, source = entry

        # Ensure strand, chrom, and source match rest for this gene
        if gene_strand is None:
            gene_strand = strand
        elif gene_strand != strand:
            die("Found gene features on more than one strand: [%s, %s]\n%s" %
                (gene_strand, strand, str(entries)))
        if gene_chrom is None:
            gene_chrom = chrom
        elif gene_chrom != chrom:
            die("Found gene features on more than one chromosome:"
                " [%s, %s]\n%s" % (gene_chrom, chrom, str(entries)))
        if gene_source is None:
            gene_source = source
        elif gene_source != source:
            print >>sys.stderr, ("Found gene features from more than one"
                                 "source: [%s, %s]\n%s" %
                                 (gene_source, source, str(entries)))

        partial_entry = (start, end)
        if component == EXON_COMPONENT:
            exons.append(partial_entry)
        elif component == CDS_COMPONENT:
            cdss.append(partial_entry)
    return gene_chrom, gene_strand, gene_source, exons, cdss

def interpret_features_as_gene(entries, min_exons=MIN_EXONS, min_cdss=MIN_CDSS):
    """Returns a set of gene-model entries based upon the given entries.

    Interprets the given set of entries in terms of an idealized gene model
    with 5' and 3' UTRs, initial, internal, and terminal exons and introns.

    :param entries: see preprocess_entries for a description
    :returns: a new set of entries based upon this gene model
              or None if the entries don't fit the model.
    """

    gene_chrom, gene_strand, gene_source, exons, cdss = \
        preprocess_entries(entries)

    # ... removing gene if without a coding exon or enough exons
    if len(exons) < min_exons or len(cdss) < min_cdss:
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
    if gene_strand == "-":
        exons.reverse()
        introns.reverse()
        cdss.reverse()

    features = []
    # Add gene features to new list, renaming components
    #   according to an idealized gene model
    # Assumes exons and cdss in sorted order (ensured by sort() above)

    # Binds local variables: gene_chrom, gene_strand, gene_source
    def make_feature(component, start, end):
        """Makes tuple feature for feature list"""
        return (gene_chrom, start, end, gene_strand, component, gene_source)

    # first exon
    features.append(make_feature(SPLICE_COMPONENTS[0], *exons[0]))
    # first intron
    if len(introns) > 0:
        features.append(make_feature(SPLICE_COMPONENTS[1], *introns[0]))
    # internal exons
    for exon in exons[1:-1]:
        features.append(make_feature(SPLICE_COMPONENTS[2], *exon))
    # internal introns
    for intron in introns[1:-1]:
        features.append(make_feature(SPLICE_COMPONENTS[3], *intron))
    # last intron
    if len(introns) > 0:
        features.append(make_feature(SPLICE_COMPONENTS[5], *introns[-1]))
    # last exon
    features.append(make_feature(SPLICE_COMPONENTS[4], *exons[-1]))


    # If there were no CDS's, then there is no point in dealing with UTR/CDS
    #   regions, so skip all that.
    if len(cdss) == 0:
        return features

    # Binds local variable gene_strand to strand-correct
    def upstream(x, y):
        """Returns true if x starts before (is 5') of y, false otherwise"""
        return (gene_strand == "+" and x[0] < y[0]) or \
            (gene_strand == "-" and x[1] > y[1])

    # Go through exons and introns and pick out those that are in UTR regions
    # Exons are a little trickier because the UTR region has to be trimmed to
    #   the CDS if they overlap
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
        if gene_strand == "+":
            if exon_start < first_cds_start:
                exon_end = min(exon_end, first_cds_start)
                utr_list = UTR5p_exons
            elif exon_end > last_cds_end:
                exon_start = max(exon_start, last_cds_end)
                utr_list = UTR3p_exons
        else:  # gene_strand == "-":
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
        features.append(make_feature(CODING_COMPONENTS[0], *UTR5p_exons[0]))
    # 5' UTR introns
    for intron in UTR5p_introns:
        features.append(make_feature(CODING_COMPONENTS[1], *intron))
    # internal 5' UTR
    for exon in UTR5p_exons[1:-1]:
        features.append(make_feature(CODING_COMPONENTS[2], *exon))
    # terminal 5' UTR
    if len(UTR5p_exons) > 0:
        features.append(make_feature(CODING_COMPONENTS[3], *UTR5p_exons[-1]))
    # first CDS
    features.append(make_feature(CODING_COMPONENTS[4], *first_cds))
    # last CDS
    features.append(make_feature(CODING_COMPONENTS[5], *last_cds))
    # initial 3' UTR
    if len(UTR3p_exons) > 0:
        features.append(make_feature(CODING_COMPONENTS[6], *UTR3p_exons[0]))
    # internal 3' UTR
    for exon in UTR3p_exons[1:-1]:
        features.append(make_feature(CODING_COMPONENTS[7], *exon))
    # 3' UTR introns
    for intron in UTR3p_introns:
        features.append(make_feature(CODING_COMPONENTS[8], *intron))
    # terminal 3' UTR
    if len(UTR3p_exons) > 0:
        features.append(make_feature(CODING_COMPONENTS[9], *UTR3p_exons[-1]))

    return features

def get_transcript_length(entries):
    """Returns the length of the transcript.

    Assumes entries are all for the same transcript

    :param entries: each entry is (chrom, start, end, strand, component, source)
                    and all entries should correspond to the same gene_id
    :type entries: list of 6-tuples
    :rtype: integer or None

    """
    min_start = None
    max_end = None
    for entry in entries:
        (chrom, start, end, strand, component, source) = entry
        if component == EXON_COMPONENT:
            if min_start is None or start < min_start:
                min_start = start
            if max_end is None or end > max_end:
                max_end = end

    if min_start is None or max_end is None:
        return None
    else:
        return max_end - min_start

def get_longest_transcript(transcript_dict):
    """Return the longest transcript in the dict.

    :type transcript_dict: transcript_id -> entries
    :returns: tuple of (transcript_id, entries) for the longest transcript
              or None if there were no transcripts.
    """
    max_length = 0
    longest = None
    for id, entries in transcript_dict.iteritems():
        length = get_transcript_length(entries)
        if length > max_length:
            max_length = length
            longest = (id, entries)
    return longest

def load_gtf_data(gtf_filename, min_exons=MIN_EXONS, min_cdss=MIN_CDSS):
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

    If an exon is a coding exon, there must also be a CDS line
    establishing the region of that exon that is translated. Assumes each CDS
    region is completely contained within a single exon.

    A gene is a coding gene iff there is at least one CDS feature in it.

    Genes that do not fit this idealized model (because they only have one
    exon or do not contain a start codon) are skipped.
    """
    # Start with establishing a dict:
    #   gene_id -> dict(transcript_id -> (chrom, start, end, strand, component))
    #   start is 0-indexed, end is non-inclusive (BED)
    gene_dict = defaultdict(partial(defaultdict, list))
    gene_id_col = 0
    transcript_id_col = 1
    with maybe_gzip_open(gtf_filename) as ifp:
        for line in ifp:
            # Ignore comment lines
            if line.startswith("#"): continue

            # Parse tokens from GTF line
            try:
                tokens = line.strip().split("\t")
                chrom = tokens[0]
                source = tokens[1]
            except ValueError:
                die("Error parsing chrom (field 1) or source (field 2)"
                    " from GTF line: %s" % line)
            try:
                start = int(tokens[3]) - 1
                end = int(tokens[4])
            except ValueError:
                die("Error parsing start (field 4) or end (field 5)"
                    " from GTF line: %s" % line)
            strand = tokens[6]
            if strand != "+" and strand != "-":
                die("Expected +/- strand, but found: %s in GTF line:"
                    " %s" % (strand, line))
            component = tokens[2]

            properties = tokens[8].split("; ")
            if not properties[gene_id_col].startswith("gene_id"):
                # Try to find the gene_id column
                gene_id_col = None
                for col_offset, prop in enumerate(properties):
                    if prop.startswith("gene_id"):
                        gene_id_col = col_offset

                if gene_id_col is None:
                    die("Could not find the gene_id column in GTF line: %s" %
                        line)
            gene_id = properties[gene_id_col].split()[1].strip('"')

            if not properties[transcript_id_col].startswith("transcript_id"):
                # Try to find the transcript_id column
                transcript_id_col = None
                for col_offset, prop in enumerate(properties):
                    if prop.startswith("transcript_id"):
                        transcript_id_col = col_offset

                if transcript_id_col is None:
                    die("Could not find the transcript_id column in"
                        " GTF line: %s" % line)
            transcript_id = properties[transcript_id_col].split()[1].strip('"')

            # Add feature to dict
            entry = (chrom, start, end, strand, component, source)
            gene_dict[gene_id][transcript_id].append(entry)

    # Eventually create feature dict: chrom -> list(features)
    data = defaultdict(list)
    for gene_id, transcript_dict in gene_dict.iteritems():

        # Select only longest transcript
        transcript_id, transcript = get_longest_transcript(transcript_dict)
        gene_model = interpret_features_as_gene(transcript,
                                                min_exons=min_exons,
                                                min_cdss=min_cdss)

        if gene_model is not None:
            # Add features to a normal, chrom-based feature dict
            for gene_part in gene_model:
#                 print str(gene_part)
                chrom, start, end, strand, component, source = gene_part
                feature = dict(start=start, end=end, strand=strand,
                               name=component, source=source)
                data[chrom].append(feature)

#         raw_input()

    # Sort features by ascending start
    for chrom_features in data.itervalues():
        chrom_features.sort(key=itemgetter("start"))

    return data


## Accepts data from dict: chr -> dict {"start", "end", "strand", "name"})
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
def calc_aggregation(segmentation, mode, features, groups, components=[],
                     component_bins=None, quick=False, by_groups=False):
    assert len(groups) > 0

    if component_bins is None:
        component_bins = get_component_bins(components=components)
    else:
        for component in components:
            assert component in component_bins

    if not by_groups:
        assert len(groups) == 1

    labels = segmentation.labels

    print >>sys.stderr, "\tGroups: %s" % groups
    print >>sys.stderr, "\tComponents and bins:"
    for component in components:
        print >>sys.stderr, "\t\t",
        # Try to substitute number of bins in first
        try:
            print >>sys.stderr, component % component_bins[component]
        except TypeError:
            print >>sys.stderr, "%s: %d" % (component,
                                            component_bins[component])

    # dict:
    #   key: feature_group
    #   value: dict:
    #      key: component_name
    #      value: numpy.array histogram [bin, label_key]
    counts = dict([(group,
                    dict([(component, fill_array(0, (bins, len(labels))))
                          for component, bins in component_bins.iteritems()]))
                   for group in groups])

    counted_features = 0

    for chrom, segments in segmentation.chromosomes.iteritems():
        print >>sys.stderr, "\t%s" % chrom

        # XXX: Segments and features are in sorted order. TAKE ADVANTAGE

        # Get bounds on segmentation for chr (since segments are sorted)
        segmentation_start = segments['start'].min()
        segmentation_end = segments['end'].max()

        # Map entire chromosome into a list of labels
        # XXXopt: THIS LINE TAKES 3.2 GB OF MEMORY!!!!
        segment_map = map_segments(segments, labels, segmentation_end)

        chrom_features = features[chrom]
        # For each feature, tally segments in window
        for feature in chrom_features:
            feature_counted = False
            if feature["end"] < segmentation_start or \
                    feature["start"] > segmentation_end:
                continue  # Feature outside segmentation

            if by_groups:
                group = get_feature_group(feature, mode)
            else:
                group = groups[0]

            group_counts = counts[group]

            # Spread internal bins throughout feature
            component_windows = calc_feature_windows(feature, mode, components,
                                                     component_bins)
            # Scan window, tallying segments observed
            for component, window in component_windows.iteritems():
                component_counts = group_counts[component]
                for bin, bp in enumerate(window):
                    try:
                        label_key = segment_map[bp]
                    except IndexError:
                        continue  # Window outside segmentation. Ignore
                    if label_key in labels:
                        component_counts[bin, label_key] += 1
                        if not feature_counted:
                            counted_features += 1
                            feature_counted = True
        if quick: break

    return (counts, counted_features)

def make_row(labels, row_data):
    values = {}
    for label_key, label in labels.iteritems():
        values[label] = row_data[label_key]

    return values

## Saves the data to a tab file
def save_tab(labels, counts, components, component_bins,
             counted_features, dirpath, clobber=False, namebase=NAMEBASE):
    metadata = {"num_features": counted_features}
    label_keys, labels = get_ordered_labels(labels)
    for label in labels.values():
        assert label not in STATIC_FIELDNAMES

    fieldnames = STATIC_FIELDNAMES + [labels[label_key]
                                      for label_key in label_keys]
    with tab_saver(dirpath, namebase, fieldnames=fieldnames, metadata=metadata,
                   clobber=clobber) as saver:
        for group in counts:
            for component in components:
                hist = counts[group][component]
                # Try to substitute component bins into component names
                try:
                    component_name = component % component_bins[component]
                except TypeError:
                    component_name = component

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
def save_plot(dirpath, num_labels, mode, spacers=len(EXON_COMPONENTS),
              clobber=False, mnemonics=[], normalize=False):
    start_R()

    tabfilename = make_tabfilename(dirpath, NAMEBASE)
    if not os.path.isfile(tabfilename):
        die("Unable to find tab file: %s" % tabfilename)

    png_size = int(PNG_BLOCK_SIZE * math.ceil(math.sqrt(num_labels)))
    # Plot data in tab file
    with image_saver(dirpath, NAMEBASE,
                     height=png_size,
                     width=png_size,
                     clobber=clobber):
        if mode != GENE_MODE:
            spacers = []

        r.plot(r["plot.aggregation"](tabfilename,
                                     mnemonics=mnemonics,
                                     genemode=(mode == GENE_MODE),
                                     spacers=spacers,
                                     normalize=normalize))

def save_html(dirpath, featurefilename, mode, num_features, groups,
              components, clobber=False, normalize=False):
    featurebasename = os.path.basename(featurefilename)
    title = "%s (%s)" % (HTML_TITLE_BASE, featurebasename)
    grouplist = list2html(groups, code=True)
    if normalize:
        yaxis = "proportion"
    else:
        yaxis = "count"

    save_html_div(HTML_TEMPLATE_FILENAME, dirpath, NAMEBASE, clobber=clobber,
                  module=MODULE, featurefilename=featurebasename,
                  numfeatures=num_features, mode=mode, grouplist=grouplist,
                  title=title, yaxis=yaxis)

def print_array(arr, tag="", type="%d"):
    fstring = "%%s:  %s, %s, ..., %s, ..., %s, %s" % tuple([type]*5)
    print >>sys.stderr, fstring % (tag,
                                   arr[0],
                                   arr[1],
                                   arr[int(arr.shape[0] / 2)],
                                   arr[-2],
                                   arr[-1])

## Package entry point
def validate(bedfilename, featurefilename, dirpath,
             flank_bins=FLANK_BINS, region_bins=REGION_BINS,
             intron_bins=INTRON_BINS, exon_bins=EXON_BINS,
             by_groups=False, mode=DEFAULT_MODE, clobber=False,
             quick=False, replot=False, noplot=False, normalize=False,
             min_exons=MIN_EXONS, min_cdss=MIN_CDSS,
             mnemonicfilename=None):
    setup_directory(dirpath)
    segmentation = load_segmentation(bedfilename)

    assert segmentation is not None

    print >>sys.stderr, "Using file %s" % featurefilename
    if mode == "gene":
        features = load_gtf_data(featurefilename,
                                 min_exons=min_exons,
                                 min_cdss=min_cdss)
    else:
        features = load_gff_data(featurefilename)
    assert features is not None

    groups = ["features"]
    # Parse components from features
    if mode == GENE_MODE:
        components = GENE_COMPONENTS
    elif mode == REGION_MODE:
        components = REGION_COMPONENTS
        if by_groups:
            groups = feature_field_values(features, "name")
    elif mode == POINT_MODE:
        components = POINT_COMPONENTS
        if by_groups:
            groups = feature_field_values(features, "name")

    labels = segmentation.labels
    mnemonics = map_mnemonics(labels, mnemonicfilename)

    num_features = 0
    for chrom_features in features.itervalues():
        num_features += len(chrom_features)
    print >>sys.stderr, "\tAggregating over %d features" % num_features

    if not replot:
        component_bins = get_component_bins(components,
                                            flank_bins=flank_bins,
                                            region_bins=region_bins,
                                            intron_bins=intron_bins,
                                            exon_bins=exon_bins)

        res = calc_aggregation(segmentation, mode=mode, features=features,
                               groups=groups, components=components,
                               component_bins=component_bins, quick=quick,
                               by_groups=by_groups)
        counts, counted_features = res
        for group in groups:
            for component in components:
                # Allow component bin substitution
                try:
                    component_name = component % component_bins[component]
                except TypeError:
                    component_name = component
                component_counts = counts[group][component].sum(axis=1)
                print_array(component_counts, tag="%s:%s count" % \
                                (group, component_name))

        save_tab(labels, counts, components, component_bins,
                 counted_features, dirpath, clobber=clobber)

    if not noplot:
        save_plot(dirpath, len(labels), mode, clobber=clobber,
                  mnemonics=mnemonics, normalize=normalize)

    save_html(dirpath, featurefilename, mode=mode, groups=groups,
              components=components, num_features=num_features,
              clobber=clobber, normalize=normalize)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTIONS] BEDFILE FEATUREFILE"
    description = ("Plot the enrichment of the BEDFILE groups relative to"
                   " the position of features in FEATUREFILE."
                   " BEDFILE is a BED4+ segmentation, where the name column"
                   " corresponds to the segment group. FEATUREFILE should"
                   " be in GFF or GTF format and can be grouped by the"
                   " feature column by supplying --groups.")
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Input options")
    group.add_option("--mnemonic-file", dest="mnemonicfilename",
                      default=None,
                      help="If specified, labels will be shown using"
                      " mnemonics found in this file")
    group.add_option("-o", "--outdir",
                      dest="outdir", default="%s" % MODULE,
                      help="File output directory (will be created"
                      " if it does not exist) [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Flags")
    group.add_option("--clobber", action="store_true",
                     dest="clobber", default=False,
                     help="Overwrite existing output files if the specified"
                     " directory already exists.")
    group.add_option("--quick", action="store_true",
                     dest="quick", default=False,
                     help="Compute values only for one chromosome.")
    group.add_option("--replot", action="store_true",
                     dest="replot", default=False,
                     help="Load data from output tab files and"
                     " regenerate plots instead of recomputing data")
    group.add_option("--noplot", action="store_true",
                     dest="noplot", default=False,
                     help="Do not generate plots")
    group.add_option("--groups", action="store_true",
                     dest="by_groups", default=False,
                     help="Separate data into different groups based upon"
                     " FEATUREFILE feature field"
                     " if --mode=region or --mode=point. Does"
                     " nothing if --mode=gene.")
    group.add_option("--normalize", action="store_true",
                     dest="normalize", default=False,
                     help="Plot the relative frequency instead of the counts"
                     " for each BEDFILE group (normalize over BEDFILE groups)")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Main aggregation options")
    group.add_option("-m", "--mode", choices=MODES,
                     dest="mode", type="choice", default=DEFAULT_MODE,
                     help="one of: "+str(MODES)+", --gene not implemented"
                     " [default: %default]")
    group.add_option("-f", "--flank-bins",
                     dest="flankbins", type="int", default=FLANK_BINS,
                     help="Aggregate this many base pairs off each"
                     " end of feature or gene [default: %default]")
    group.add_option("-r", "--region-bins", type="int",
                     dest="regionbins", default=REGION_BINS,
                     help="If --mode=region, aggregate over each internal"
                     "feature using this many evenly-spaced bins"
                     " [default: %default]")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Gene aggregation options")
    group.add_option("-i", "--intron-bins", type="int",
                     dest="intronbins", default=INTRON_BINS,
                     help="Aggregate over each intron"
                     "using this many evenly-spaced bins"
                     " [default: %default]")
    group.add_option("-e", "--exon-bins", type="int",
                     dest="exonbins", default=EXON_BINS,
                     help="Aggregate over each exon"
                     "using this many evenly-spaced bins"
                     " [default: %default]")
#     group.add_option("--min-exons", type="int",
#                      dest="min_exons", default=MIN_EXONS,
#                      help="Only consider genes with at least this many exons"
#                      " [default: %default]")
#     group.add_option("--min-coding-exons", type="int",
#                      dest="min_cdss", default=MIN_CDSS,
#                      help="Only consider genes with at least this many coding"
#                      " exon [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(args)

    if len(args) < 2:
        parser.error("Insufficient number of arguments")

    if options.noplot and options.replot:
        parser.error("noplot and replot are contradictory")

    return (options, args)

## Command-line entry point
def main(args=sys.argv[1:]):
    (options, args) = parse_options(args)
    bedfilename = args[0]
    featurefilename = args[1]
    kwargs = {"flank_bins": options.flankbins,
              "region_bins": options.regionbins,
              "intron_bins": options.intronbins,
              "exon_bins": options.exonbins,
              "clobber": options.clobber,
              "quick": options.quick,
              "replot": options.replot,
              "noplot": options.noplot,
              "by_groups": options.by_groups,
              "normalize": options.normalize,
              "mode": options.mode,
#               "min_exons": options.min_exons,
#               "min_cdss": options.min_cdss,
              "mnemonicfilename": options.mnemonicfilename}
    validate(bedfilename, featurefilename, options.outdir, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
