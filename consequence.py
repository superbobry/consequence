# -*- coding: utf-8 -*-
"""
     consequence
     ~~~~~~~~~~~

     Find the consequences of a given genomic sequence; no pun intended.
"""

import cPickle
import operator
import os.path
import sys
from collections import defaultdict, namedtuple

import vcf
import pysam


#: Path to global VCF index.
CONSEQUENCE_INDEX = "consequence.index"

_Chunk = namedtuple("_Chunk", list("ACGT"))

class Chunk(_Chunk):
    def __new__(cls, A=None, C=None, G=None, T=None):
        return super(cls, Chunk).__new__(cls, A, C, G, T)


def update_index(seq_path, seq_id):
    if os.path.exists(CONSEQUENCE_INDEX):
        index = cPickle.load(open(CONSEQUENCE_INDEX, "rb"))
    else:
        # Yeah-yeah, it's not a 'persistent trie', because there isn't
        # one avaiable for Python yet.
        index = {}

    for record in vcf.VCFReader(open(seq_path, "rb")):
        if not record.is_snp:
            continue

        index.setdefault(record.POS, Chunk(**{record.REF: None}))

        for alt in filter(operator.truth, record.ALT):
            known = getattr(index[record.POS], alt) or set()

            if seq_id not in known:
                known.add(seq_id)
                index[record.POS] = index[record.POS]._replace(**{alt: known})

    cPickle.dump(index, open(CONSEQUENCE_INDEX, "wb"))


def naive_lookup(seq_path):
    index = cPickle.load(open(CONSEQUENCE_INDEX, "rb"))

    # A mapping of (pos, base) -> frequency; would be nice to store
    # read quality as well for further SNP evaluation.
    snps = defaultdict(int)
    for record in pysam.Samfile(seq_path, "rb"):
        if record.is_duplicate or record.is_unmapped:
            continue

        for pos in record.positions:
            base = record.query[pos - record.pos - record.qstart]

            if pos not in index or base not in "ACGT":
                continue

            if getattr(index[pos], base) is not None:
                snps[pos, base] += 1

    # Filter the resulting mapping and output a set of the corresponding
    # genome identifiers.
    return snps


if __name__ == "__main__":
    # update_index(seq_path, seq_id)
    # naive_lookup(seq_path)
    pass
