#!/usr/bin/env python
from __future__ import division

# Copyright 2008 Michael M. Hoffman <mmh1@washington.edu>

import sys

FIELDNAMES = ["seqname", "source", "feature", "start", "end",  # required
              "score", "strand", "frame", "group"]

class Datum(object):
    def __init__(self, words, *args, **kwargs):
        self.__dict__ = dict(zip(FIELDNAMES, words))
        self._words = tuple(words)

    def __repr__(self):
        return "%s%s" % (self.__class__.__name__, self._words)

class NativeDatum(Datum):
    def __init__(self, *args, **kwargs):
        Datum.__init__(self, *args, **kwargs)

        # convert to zero-based, half-open:
        # http://genome.ucsc.edu/FAQ/FAQformat#format3
        try:
            self.start = int(self.start) - 1
            self.end = int(self.end)
            try:
                self.score = float(self.score)
            except AttributeError:
                pass
            except ValueError:
                pass

            # attributes: http://genome.ucsc.edu/FAQ/FAQformat#format4
            if "gtf" in kwargs and kwargs["gtf"]:
                self.attributes = {}
                for attribute in self.group.rstrip("; ").split("; "):
                    type, value = attribute.split(" ", 1)
                    self.attributes[type] = value.strip('"')  # Remove quotes

                assert "gene_id" in self.attributes

                # New gencode annotations rows (>=v24) with featuretype 'gene'
                # no longer have a 'transcript_id' attribute. 
                # Add the 'transcript_id' key to the attribute dict and set its
                # value to the the value of 'gene_id' like it was in older
                # gencode annotations.
                if self.feature == "gene" and \
                   "transcript_id" not in self.attributes:
                   
                    self.attributes["transcript_id"] = self.attributes["gene_id"]

                assert "transcript_id" in self.attributes                 

            self._words = ((self.seqname, self.start, self.end)
                           + self._words[3:])
        except:
            print >>sys.stderr, "Error processing line: %s" % \
                  "\t".join(self._words)
            raise

def read(iterator, datum_cls=Datum, *args, **kwargs):
    try:
        for line in iterator:
            line = line.split("#", 1)[0].rstrip() # Ignore comments
            if not line: continue

            words = line.split("\t")  # Tab-delimited
            if words[0] == "track":  # Ignore any track lines
                continue

            assert len(words) >= 5

            yield datum_cls(words, *args, **kwargs)
    except:
        print >>sys.stderr, "Error processing line: %s" % line
        raise

def read_native(*args, **kwargs):
    return read(datum_cls=NativeDatum, *args, **kwargs)

def main(args=sys.argv[1:]):
    pass

if __name__ == "__main__":
    sys.exit(main())
