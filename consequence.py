# -*- coding: utf-8 -*-
"""
     consequence
     ~~~~~~~~~~~

     Find the consequences of a given genomic sequence; no pun intended.
"""

from __future__ import print_function

import cPickle
import glob
import operator
import os.path
import sys
from collections import defaultdict, namedtuple, MutableMapping

import vcf
import pysam


_Chunk = namedtuple("_Chunk", list("ACGT"))

class Chunk(_Chunk):
    def __new__(cls, A=None, C=None, G=None, T=None):
        return super(cls, Chunk).__new__(cls, A, C, G, T)


class Genome(MutableMapping):
    def __init__(self, base_path=None):
        if base_path and not os.path.isdir(base_path):
            raise RuntimeError("Invalid genome base path: {0!r}"
                               .format(base_path))

        self.base_path = base_path or os.getcwd()
        self.cache = {}

    def dump(self):
        for chr in self.cache:
            self.dump_index(chr)

    def load(self):
        for chr in glob.glob(os.path.join(self.base_path, "*.index")):
            self.load_index(chr)

    def load_index(self, chr):
        path = os.path.join(self.base_path, "{0}.index".format(chr))
        if os.path.isfile(path):
            # Note(Sergei): this is probably the dumbest thing possible,
            # but hey -- it's still better than nothing ;)
            self.cache[chr] = cPickle.load(open(path, "rb"))
        else:
            # Yeah-yeah, it's not a 'persistent trie', because there
            # isn't one avaiable for Python yet.
            self.cache[chr] = {}

    def dump_index(self, chr):
        path = os.path.join(self.base_path, "{0}.index".format(chr))
        cPickle.dump(self.cache[chr], open(path, "wb"))

    def __len__(self):
        # Obviously, this does *NOT* include non-cached indices.
        return sum(map(len, self.cache))

    def __iter__(self):
        return ((chr, pos) for pos in self.cache[chr] for chr in self.cache)

    def __contains__(self, (chr, pos)):
        if chr not in self.cache:
            self.load_index(chr)

        return pos in self.cache[chr]

    def __getitem__(self, (chr, pos)):
        if chr not in self.cache:
            self.load_index(chr)

        return self.cache[chr][pos]

    def __setitem__(self, (chr, pos), chunk):
        if chr not in self.cache:
            self.load_index(chr)

        self.cache[chr][pos] = chunk

    def __delitem__(self, (chr, pos)):
        del self.cache[chr][pos]


def update_index(seq_path, seq_id):
    g = Genome()

    for record in vcf.VCFReader(open(seq_path, "rb")):
        if not record.is_snp:
            continue

        g.setdefault((record.CHROM, record.POS), Chunk(**{record.REF: None}))

        for alt in filter(operator.truth, record.ALT):
            known = getattr(g[record.CHROM, record.POS], alt) or set()

            if seq_id not in known:
                known.add(seq_id)
                g[record.CHROM, record.POS] = \
                    g[record.CHROM, record.POS]._replace(**{alt: known})

    g.dump()


def naive_lookup(seq_path):
    g = Genome()
    f = pysam.Samfile(seq_path, "rb")

    # A mapping of (pos, base) -> frequency; would be nice to store
    # read quality as well for further SNP evaluation.
    snps = defaultdict(int)
    for record in f:
        if record.is_duplicate or record.is_unmapped or record.tid < 0:
            continue

        chrom = f.getrname(record.tid)

        for pos in record.positions:
            base = record.query[pos - record.pos - record.qstart]

            if base not in "ACGT" or (chrom, pos) not in g:
                continue

            if getattr(g[chrom, pos], base) is not None:
                snps[pos, base] += 1

    # Filter the resulting mapping and output a set of the corresponding
    # genome identifiers.

    # Should we store frequency of an SNP in a particular sample in the
    # index?
    return snps


if __name__ == "__main__":
    # update_index(seq_path, seq_id)
    # naive_lookup(seq_path)
    pass
