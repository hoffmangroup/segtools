#!/usr/bin/env python
from __future__ import division, with_statement

"""
html.py

HTML utilities for segtools
"""

import os
import re
import sys
import time

from functools import partial
from pkg_resources import resource_string
from shutil import copy
from string import Template

from . import log, Segmentation, die, add_common_options
from .common import check_clobber, get_ordered_labels, make_divfilename, \
     make_id, make_filename, make_tabfilename, map_mnemonics, NICE_EXTS, \
     PKG_RESOURCE
from .version import __version__

MNEMONIC_TEMPLATE_FILENAME = "mnemonic_div.tmpl"
HEADER_TEMPLATE_FILENAME = "html_header.tmpl"
FOOTER_TEMPLATE_FILENAME = "html_footer.tmpl"
GENOMEBROWSER_URL = "http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=track"
GENOMEBROWSER_OPTIONS = {"autoScale":"off",
                         "viewLimits":"0:1",
                         "visibility":"full",
                         "itemRgb":"on",
                         "name":"segtools"}
GENOMEBROWSER_LINK_TMPL = """
<li>Link to view this segmentation in the UCSC genome browser:<br />
<script type="text/javascript">print_genomebrowser_link("%s");</script>
</li>"""

DESCRIPTION_MODULE = ("description", "Segmentation information")
MNEMONIC_MODULE = ("mnemonics", "Mnemonics")

template_string = partial(resource_string, PKG_RESOURCE)

def template_substitute(filename):
    """
    Simplify import resource strings in the package
    """
    return Template(template_string(filename)).safe_substitute

def tuple2link(entry):
    """entry should be a (url, text) tuple"""
    return '<a href="%s">%s</a>' % entry

def list2html(list, code=False, link=False):
    """
    If link is True: each element of the list should be a (divID, label) tuple
    """
    result = ["<ul>"]
    entrystr = "%s"
    if link:
        entrystr = '<a href="#%s">%s</a>'
    if code:
        entrystr = "<code>%s</code>"

    entrystr = "<li>%s</li>" % entrystr
    for entry in list:
        result.append(entrystr % entry)

    result.append("</ul>")
    return "\n".join(result)

def tab2html(tabfile, header=True, mnemonicfile=None):
    """
    Given a tab file table with a header row, generates an html string which
    for pretty-display of the table data.
    """
    if not os.path.isfile(tabfile):
        return "<File not found>"

    result = []
    result.append('\n<table border="1" cellpadding="4" cellspacing="1">')
    # if the tabfile exists, write in htmlhandle an html table
    with open(tabfile) as ifp:
        # Read past comments
        line = ifp.readline()
        while line.startswith("#"):
            line = ifp.readline()

        if header:
            # Write colname row from header row of file
            fields = line.split("\t")
            result.append("<tr>")

            for f in fields:
                result.append('<td style="background-color:'
                             'rgb(204, 204, 204)">%s</td>' % f)

            result.append("</tr>")
            line = ifp.readline()

        lines = [line] + ifp.readlines()
        rows = [line.split("\t") for line in lines]
        row_names = [row[0] for row in rows]
        # Make basic labels for row names
        row_order = range(0, len(rows))
        row_labels = dict(zip(row_order, row_names))
        if mnemonicfile:
            # Substitute these labels with mnemonics
            mnemonics = map_mnemonics(row_labels, mnemonicfile)
            row_order, row_labels = get_ordered_labels(row_labels,
                                                       mnemonics)

        for row_key in row_order:
            entry = ["<tr>"]
            fields = rows[row_key]
            fields[0] = row_labels[row_key]  # Replace row_name
            for f in fields:
                entry.append("<td>%s</td>" % f)

            entry.append("</tr>")
            result.append("".join(entry))

    result.append("</table>\n")
    return "\n".join(result)

def find_output_files(dirpath, namebase, d={}, tag=""):
    exts = NICE_EXTS
    # Add filenames of common present files to dict
    for extname, ext in exts.iteritems():
        filename = make_filename(dirpath, namebase, ext)
        if os.path.isfile(filename):
            key = "%s%sfilename" % (tag, extname)
            assert key not in d
            d[key] = filename

    return d

def form_template_dict(dirpath, namebase, module=None, extra_namebases={},
                       mnemonicfile=None, tables={}, **kwargs):
    """
    Given information about the current validation, generates a dictionary
    suitable for HTML template substitution.

    The output directory (dirpath) is searched for files of the form:
    <namebase>.<ext> for common exts. If found, the filename is linked
    in under the variable <ext>filename

    extra_namebases: a dict of tag -> namebase string
    For each extra namebase, any found files will be linked under
    <tag><ext>filename, as opposed to the main namebase, which is just
    under <ext>filename.

    module: The name of the module generating the dict. If specified,
    an id variable will be generated based upon the module and dirpath
    (a pseudo-unique identifier for the div file).

    tables: a dict from tag -> to table file namebase string.
    The rownames of the tables will try to be substituted with the mnemonics
    in mnemonicfile, if provided. Variable of the form <tag>table and
    <tag>tablefilename will be created with the table HTML and table file
    name, respectively. The dict value can also be a tuple of two elements,
    where the first is the namebase string and the second is a string
    mode. Currently, only None and "exact" are supported,
    with the latter causing the table to not be mnemonic-substituted.

    Any other keyword args supplied are linked into the dictionary.
    """
    # Find default files for all namebases
    d = {}
    find_output_files(dirpath, namebase, d=d)
    for tag, nb in extra_namebases.iteritems():
        find_output_files(dirpath, nb, d=d, tag=tag)

    if module is not None:
        arg = "id"
        assert arg not in d
        d[arg] = make_id(module, dirpath)

    # Add any tables
    for tag, table in tables.iteritems():
        if isinstance(table, tuple):
            table, tablemode = table
        else:
            tablemode = None

        tablefilename = make_tabfilename(dirpath, table)
        filearg = "%stablefilename" % tag
        tablearg = "%stable" % tag
        if tablemode is None:
            val = tab2html(tablefilename, mnemonicfile=mnemonicfile)
        elif tablemode == "exact":
            val = tab2html(tablefilename)
        else:
            raise ValueError("Invalid table mode: %s" % tablemode)

        assert filearg not in d and tablearg not in d  # Don't overwrite
        d[filearg] = tablefilename
        d[tablearg] = val

    for arg, val in kwargs.iteritems():
        assert arg not in d  # Don't overwrite
        d[arg] = val

    return d

def write_html_div(dirpath, namebase, html, clobber=False):
    """
    Write the given html div string to an appropriate file
    """
    filename = make_divfilename(dirpath, namebase)
    check_clobber(filename, clobber)

    with open(filename, "w") as ofp:
        ofp.write("%s\n" % html)

def save_html_div(template_filename, dirpath, namebase,
                  clobber=False, verbose=True, **kwargs):
    """Save an HTML div file for a module run by subsituting a template file"""
    fields = form_template_dict(dirpath, namebase, **kwargs)
    div_filename = make_divfilename(dirpath, namebase)
    log("Saving html data to file: %s" % div_filename, verbose)
    try:
        html = template_substitute(template_filename)(fields)
    except KeyError, e:
        log("Error: Missing data: %s. Skipping HTML output." % e)
        return

    write_html_div(dirpath, namebase, html, clobber=clobber)

def form_mnemonic_div(mnemonicfile, results_dir, clobber=False, verbose=True):
    """Copy mnemonic file to results_dir and create the HTML div with a link"""
    filebase = os.path.basename(mnemonicfile)
    link_file = os.path.join(results_dir, filebase)

    if os.path.exists(link_file) and os.path.samefile(mnemonicfile, link_file):
        log("Mnemonic file already in place: %s" % mnemonicfile)
    else:
        check_clobber(link_file, clobber)

        try:
            copy(mnemonicfile, link_file)
        except (IOError, os.error):
            log("Error: could not copy %s to %s. Linking"
                " the the former." % (mnemonicfile, link_file))
            link_file = mnemonicfile  # Link to the existing file
        else:
            log("Copied %s to %s" % (mnemonicfile, link_file), verbose)

    fields = {}
    fields["tabfilename"] = link_file
    fields["table"] = tab2html(mnemonicfile)
    div = template_substitute(MNEMONIC_TEMPLATE_FILENAME)(fields)
    return div

def make_genomebrowser_url(options, urltype):
    """Makes URL for genomebrowser (minus javascript-added file path)

    urltype: either "data" or "bigData"

    """
    url = GENOMEBROWSER_URL
    for k, v in options.iteritems():
        url += " %s=%s " % (k, v)
    url += " %sUrl=" % urltype
    return url

def form_html_header(bedfilename, modules, layeredbed=None, bigbed=None):
    segtool, segtracks = Segmentation.get_bed_metadata(bedfilename)
    bedfilebase = os.path.basename(bedfilename)
    fields = {}
    fields["bedfilename"] = bedfilebase
    fields["numsegtracks"] = len(segtracks)
    fields["segtracks"] = list2html(segtracks, code=True)
    fields["segtool"] = segtool
    fields["bedmtime"] = time.strftime(
        "%m/%d/%Y %I:%M:%S %p", time.localtime(os.path.getmtime(bedfilename)))
    fields["modules"] = list2html(modules, link=True)
    fields["otherbeds"] = ""
    fields["genomebrowserurl"] = ""
    fields["genomebrowserlink"] = ""

    otherbeds = []
    if layeredbed:
        try:
            layeredbed = os.path.relpath(layeredbed)
        except AttributeError:
            pass
        otherbeds.append(tuple2link((layeredbed, "layered")))
    if bigbed:
        try:
            bigbed = os.path.relpath(bigbed)
        except AttributeError:
            pass
        otherbeds.append(tuple2link((bigbed, "bigBed")))

    if layeredbed or bigbed:
        options = GENOMEBROWSER_OPTIONS
        options["description"] = bedfilebase
        # Specify type (only) if using bigBed
        if bigbed:
            options["type"] = "bigBed"
            urltype = "bigData"
            datafile = bigbed
        else:
            urltype = "data"
            datafile = layeredbed
            if "type" in options:
                del options["type"]

        # Specify genomebrowser values to substitute
        fields["genomebrowserurl"] = make_genomebrowser_url(options, urltype)
        fields["genomebrowserlink"] = GENOMEBROWSER_LINK_TMPL % datafile

    if len(otherbeds) > 0:
        fields["otherbeds"] = "(%s)" % (", ".join(otherbeds))

    header = template_substitute(HEADER_TEMPLATE_FILENAME)(fields)
    return header

def form_html_footer():
    return template_substitute(FOOTER_TEMPLATE_FILENAME)()

def find_divs(rootdir=os.getcwd(), verbose=True):
    """Look one level deep in directory, adding any .div files found.
    """
    divs = []
    log("Searching directories under %s for files ending in .div..." % \
        rootdir)
    for foldername in sorted(os.listdir(rootdir)):
        folderpath = os.path.join(rootdir, foldername)
        if os.path.isdir(folderpath):
            for filename in sorted(os.listdir(folderpath)):
                root, ext = os.path.splitext(filename)
                if ext == ".div":
                    filepath = os.path.join(folderpath, filename)
                    divs.append(filepath)
                    log("Found div: %s" % filename, verbose)

    if not divs:
        log("Found no files matching %s/*/*.div" % rootdir)

    return divs

def make_html_report(bedfilename, results_dir, outfile, mnemonicfile=None,
                     clobber=False, verbose=True,
                     layeredbed=None, bigbed=None):
    check_clobber(outfile, clobber)

    divs = find_divs(results_dir, verbose=verbose)
    if len(divs) == 0:
        die("Make sure to run this from the parent directory of the"
            " module output directories or specify the --results-dir option")

    body = []
    modules = [DESCRIPTION_MODULE]
    if mnemonicfile is not None:
        modules.append(MNEMONIC_MODULE)
        div = form_mnemonic_div(mnemonicfile, results_dir, clobber=clobber,
                                verbose=verbose)
        body.append(div)

    regex = re.compile('"module" id="(.*?)".*?<h.>.*?</a>\s*(.*?)\s*</h.>',
                       re.DOTALL)
    for div in divs:
        with open(div) as ifp:
            divstring = "".join(ifp.readlines())
            matching = regex.search(divstring)
            assert matching
            module = (matching.group(1), matching.group(2))
            modules.append(module)
            body.append(divstring)

    header = form_html_header(bedfilename, modules,
                              layeredbed=layeredbed,
                              bigbed=bigbed)
    footer = form_html_footer()

    components = [header] + body + [footer]
    separator = "<br /><hr>"
    with open(outfile, "w") as ofp:
        ofp.write("%s\n" % separator.join(components))

def parse_args(args):
    from optparse import OptionGroup, OptionParser
    usage = "%prog [OPTIONS] SEGMENTATION"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)

    group = OptionGroup(parser, "Flags")
    add_common_options(group, ['clobber', 'quiet'])
    parser.add_option_group(group)

    group = OptionGroup(parser, "Linking")
    group.add_option("-m", "--mnemonic-file", dest="mnemonicfile",
                     default=None, metavar="FILE",
                     help="If specified, this mnemonic mapping will be"
                     " included in the report (this should be the same"
                     " mnemonic file used by the individual modules).")
    group.add_option("-L", "--layered-bed", dest="layeredbed",
                     default=None, metavar="FILE",
                     help="If specified, this layered BED file will be linked"
                     " into the the HTML document (assumed to be the same"
                     " data as in SEGMENTATION)")
    group.add_option("-B", "--big-bed", dest="bigbed",
                     default=None, metavar="FILE",
                     help="If specified, this bigBed file will be linked into"
                     " the the HTML document and a UCSC genome brower link"
                     " will be generated for it (assumed to be the same data"
                     " as in SEGMENTATION)")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Output")
    group.add_option("--results-dir", dest="resultsdir",
                     default=".", metavar="DIR",
                     help="This should be the directory containing all the"
                     " module output directories (`ls` should return things"
                     " like \"length_distribution/\", etc)"
                     " [default: %default]")
    group.add_option("-o", "--outfile", metavar="FILE",
                     dest="outfile", default="index.html",
                     help="HTML report file (must be in current directory"
                     " or the links will break [default: %default]")
    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if options.outfile != os.path.basename(options.outfile):
        parser.error("Output file must be in current directory"
                     " (otherwise the reource paths get all messed up)")

    if len(args) < 1:
        parser.error("Insufficient number of arguments")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    bedfilename = args[0]
    kwargs = dict(options.__dict__)

    outfile = kwargs.pop('outfile')
    results_dir = kwargs.pop('resultsdir')
    make_html_report(bedfilename, results_dir, outfile, **kwargs)

if __name__ == "__main__":
    sys.exit(main())
