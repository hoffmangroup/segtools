#!/usr/bin/env python
from __future__ import division, with_statement

__version__ = "$Revision$"

"""
html.py

HTML utilities for segtools
"""

import os
import re
import sys
import time

from .common import check_clobber, die, get_bed_metadata, make_divfilename, make_id, make_filename, make_tabfilename, NICE_EXTS, template_substitute

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
def tuple2link(entry):
    """entry should be a (url, text) tuple"""
    return '<a href="%s">%s</a>' % entry

def list2html(list, code=False, link=False):
    """
    If link is True: each element of the list should be a (divID, label) tuple
    """
    result = "<ul>"
    entrystr = "%s"
    if link:
        entrystr = '<a href="#%s">%s</a>'
    if code:
        entrystr = "<code>%s</code>"

    entrystr = "<li>%s</li>" % entrystr
    for entry in list:
        result += entrystr % entry
    return result + "</ul>"

def tab2html(tabfile, header=True):
    """
    Given a tab file table with a header row, generates an html string which
    for pretty-display of the table data.
    """
    result = '\n<table border="1" cellpadding="4" cellspacing="1">'
    # if the tabfile exists, write in htmlhandle an html table
    with open(tabfile) as ifp:
        if header:
            # first read the first line, and write the header row
            line = ifp.readline()
            fields = line.split("\t")
            result += "<tr>"
            for f in fields:
                result += '<td style="background-color: \
rgb(204, 204, 204)">%s</td>' % f
            result += "</tr>"

        for line in ifp:
            result += "<tr>"
            fields = line.split("\t")
            for f in fields:
                result += "<td>%s</td>" % f
            result += "</tr>"
    result += "</table>\n"
    return result

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

def form_template_dict(dirpath, namebase, module=None,
                       extra_namebases={}, **kwargs):
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

    Any keyword arguments of the form <prefix>tablenamebase
    will be assumed to refer to a table tab files and will be converted to
    html form and linked in under a variable named <prefix>table.
    The corresponding tab file is linked under a variable named
    <prefix>tablefilename.

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

    for arg, val in kwargs.iteritems():
        index = arg.rfind("tablenamebase")
        if index >= 0: # arg ends with "tablenamebase"
            tablefilename = make_tabfilename(dirpath, val)
            if os.path.isfile(tablefilename):
                # Save html table under variable
                prefix = arg[0:index]
                arg = "%stable" % prefix
                val = tab2html(tablefilename)
                # Save table tab file under separate variable
                filearg = "%stablefilename" % prefix
                assert filearg not in d
                d[filearg] = tablefilename
            # else leave arg and val the way they were
        assert arg not in d  # Make sure wer're not overwriting anything
        d[arg] = val

    return d

def write_html_div(dirpath, namebase, html, clobber=False):
    """
    Write the given html div string to an appropriate file
    """
    filename = make_divfilename(dirpath, namebase)
    check_clobber(filename, clobber)

    with open(filename, "w") as ofp:
        print >>ofp, html

def save_html_div(template_filename, dirpath, namebase,
                  clobber=False, **kwargs):
    fields = form_template_dict(dirpath, namebase, **kwargs)

    try:
        html = template_substitute(template_filename)(fields)
    except KeyError, e:
        print >>sys.stderr, "Skipping HTML output. Missing data: %s" % e
        return

    write_html_div(dirpath, namebase, html, clobber=clobber)

def form_mnemonic_div(mnemonicfile):
    fields = {}
    fields["tabfilename"] = mnemonicfile
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
    segtool, segtracks = get_bed_metadata(bedfilename)
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

def find_divs(rootdir=os.getcwd()):
    """Look one level deep in directory, adding any .div files found.
    """
    divs = []
    for foldername in sorted(os.listdir(rootdir)):
        folderpath = os.path.join(rootdir, foldername)
        if os.path.isdir(folderpath):
            for filename in sorted(os.listdir(folderpath)):
                root, ext = os.path.splitext(filename)
                if ext == ".div":
                    filepath = os.path.join(folderpath, filename)
                    divs.append(filepath)
                    print >>sys.stderr, "Found div: %s" % filename
    return divs

def parse_args(args):
    from optparse import OptionParser
    usage = "%prog [OPTIONS] BEDFILE"
    version = "%%prog %s" % __version__
    parser = OptionParser(usage=usage, version=version)


    parser.add_option("--clobber", action="store_true",
                      dest="clobber", default=False,
                      help="Overwrite existing output files if the specified"
                      " directory already exists.")
    parser.add_option("--mnemonic-file", dest="mnemonicfile", default=None,
                      help="If specified, this mnemonic mapping will be"
                      " included in the report (this should be the same"
                      " mnemonic file used by the individual modules)")
    parser.add_option("-L", "--layered-bed", dest="layeredbed", default=None,
                      help="If specified, this layered BED file will be linked"
                      " into the the HTML document (assumed to be the same"
                      " data as in BEDFILE)")
    parser.add_option("-B", "--big-bed", dest="bigbed", default=None,
                      help="If specified, this bigBed file will be linked into"
                      " the the HTML document and a UCSC genome brower link"
                      " will be generated for it (assumed to be the same data"
                      " as in BEDFILE)")
    parser.add_option("--results-dir", dest="resultsdir", default=".",
                      help="This should be the directory containing all the"
                      " module output directories (`ls` should return things"
                      " like \"length_distribution/\", etc)"
                      " [default: %default]")
    parser.add_option("-o", "--outfile",
                      dest="outfile", default="index.html",
                      help="HTML report file (must be in current directory"
                      " or the links will break [default: %default]")

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
    mnemonicfile = options.mnemonicfile
    outfile = options.outfile
    check_clobber(outfile, options.clobber)

    divs = find_divs(options.resultsdir)
    if len(divs) == 0:
        die("Unable to find any module .div files."
            " Make sure to run this from the parent directory of the"
            " module output directories or specify the --results-dir option")

    body = []
    modules = []
    if mnemonicfile is not None:
        body.append(form_mnemonic_div(mnemonicfile))

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
                              layeredbed=options.layeredbed,
                              bigbed=options.bigbed)
    footer = form_html_footer()

    components = [header] + body + [footer]
    separator = "<br /><hr>"
    with open(outfile, "w") as ofp:
        print >>ofp, separator.join(components)

if __name__ == "__main__":
    sys.exit(main())
