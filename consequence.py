#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
     consequence
     ~~~~~~~~~~~

     Index-backed SNP caller.
"""

from __future__ import print_function

import csv
import cPickle
import glob
import itertools
import operator
import os.path
import sys
from collections import defaultdict, namedtuple, MutableMapping

import numpy as np
import opster
import pysam
import vcf
from Bio import SeqIO


_Chunk = namedtuple("_Chunk", list(["id", "A", "C", "G", "T"]))

class Chunk(_Chunk):
    def __new__(cls, id=None, A=None, C=None, G=None, T=None):
        return super(cls, Chunk).__new__(cls, id, A, C, G, T)


class Genome(MutableMapping):
    def __init__(self, base_path=None):
        if base_path and not os.path.isdir(base_path):
            raise RuntimeError("Invalid genome base path: {0!r}"
                               .format(base_path))

        self.base_path = base_path or os.getcwd()
        self.cache = {}

    def dump(self):
        for chrom in self.cache:
            path = os.path.join(self.base_path, "{0}.index".format(chrom))
            cPickle.dump(self.cache[chrom], open(path, "wb"))

    def load(self, *chroms):
        if not chroms:
            fnames = glob.glob(os.path.join(self.base_path, "*.index"))
            chroms = [os.path.basename(name)
                      for name, _ext in map(os.path.splitext, fnames)]

        for chrom in chroms:
            if chrom in self.cache:
                continue  # Oops, already got that sequence loaded.

            path = os.path.join(self.base_path, "{0}.index".format(chrom))
            if os.path.isfile(path):
                # Note(Sergei): this is probably the dumbest thing possible,
                # but hey -- it's still better than nothing ;)
                self.cache[chrom] = cPickle.load(open(path, "rb"))
            else:
                # Yeah-yeah, it's not a 'persistent trie', because there
                # isn't one avaiable for Python yet.
                self.cache[chrom] = {}

        return [self.cache[chrom] for chrom in chroms]

    def __repr__(self):
        return "<Genome: {0!r}>".format(self.cache)

    def __len__(self):
        # Obviously, this does *NOT* include non-cached indices.
        return sum(map(len, self.cache))

    def __iter__(self):
        return ((chrom, pos)
                for chrom in self.cache for pos in self.cache[chrom])

    def __contains__(self, (chrom, pos)):
        if chrom not in self.cache:
            self.load(chrom)

        return pos in self.cache[chrom]

    def __getitem__(self, (chrom, pos)):
        if chrom not in self.cache:
            self.load(chrom)

        return self.cache[chrom][pos]

    def __setitem__(self, (chrom, pos), chunk):
        if chrom not in self.cache:
            self.load(chrom)

        self.cache[chrom][pos] = chunk

    def __delitem__(self, (chrom , pos)):
        del self.cache[chrom][pos]


options = [("", "dbsnp", None, "path to dbSNP release in VCF format")]

@opster.command(options, name="index", usage="path/to/snps.vcf id")
def update_index(seq_path, seq_id, dbsnp=None, index_root=None, quiet=None):
    """Update 'consequence' index with SNPs from a given VCF file.

    If `dbsnp` argument is provided, SNPs missing in `dbsnp` will *not*
    be indexed.
    """
    g = Genome(base_path=index_root)

    dbsnp = {} if not dbsnp else \
        dict((r.POS, r.ID) for r in vcf.VCFReader(open(dbsnp, "rb")))

    for record in vcf.VCFReader(open(seq_path, "rb")):
        if not record.is_snp:
            continue

        # XXX here comes the punchline, VCF uses 1-based indexing,
        # so does SAM / BAM, but 'pysam' converts *all* indexes to
        # 0-based!
        chrom, pos = record.CHROM, record.POS - 1

        if dbsnp and pos not in dbsnp:
            continue

        g.setdefault((chrom, pos),  # vvv -- dnSNP name defaults to '*'.
                     Chunk(dbsnp.get(pos, "*"), **{record.REF: False}))

        for alt in filter(operator.truth, record.ALT):
            known = getattr(g[chrom, pos], alt) or set()
            if seq_id in known:
                continue

            known.add(seq_id)
            g[chrom, pos] = g[chrom, pos]._replace(**{alt: known})

    g.dump()


options = [("", "insert-size", 1024, "maximum expected insert size")]

@opster.command(options, name="partial", usage="path/to/reference.fasta")
def build_partial_reference(ref_path, insert_size=1024,
                            index_root=None, quiet=None):
    """Build a partial reference from a full reference sequence."""
    if not os.path.isfile(ref_path):
        raise RuntimeError("Invalid reference sequence: {0!r}"
                           .format(ref_path))

    g = Genome(base_path=index_root)
    out = open("{0}.cropped{1}".format(*os.path.splitext(ref_path)), "w")

    # Note(Sergei): it would be nice to process ref. seq. in chunks,
    # but unfortunately 'SeqIO' doesn't allow that.
    for chrom in SeqIO.parse(open(ref_path), "fasta"):
        chunks = []         # A list of chromosome chunks to save.
        left, right = 0, 0  # Crop region bounds.

        variants = sorted(*g.load(chrom.id), reverse=True)
        while variants:
            pos = variants.pop()

            if pos - right > insert_size or not variants:
                chunks.append((max(0, left - insert_size),
                               min(right + insert_size, len(chrom))))
                right = left = pos
            else:
                right = pos

        for left, right in chunks:
            # The original genomic region is encoded in sequence id.
            seq = chrom[left:right]
            seq.id = "{0}|{1}:{2}".format(seq.id, left, right)
            SeqIO.write(seq, out, "fasta")


options = [("d", "diploid", None, "output results for a diploid genome")]

@opster.command(options, name="lookup", usage="path/to/aligned_reads.bam")
def naive_lookup(seq_path, diploid=False, index_root=None, quiet=None):
    """Lookup SNPs from aligned reads in the 'consequnce' index."""
    g = Genome(base_path=index_root)

    # SAM / BAM file *must* be aligned to *partial* reference, generated
    # by 'build_partial_reference'.
    f = pysam.Samfile(seq_path, "rb")

    snps = defaultdict(lambda: defaultdict(int))
    for record in f:
        if record.is_duplicate or record.is_unmapped or record.tid < 0:
            continue

        chrom, mark = f.getrname(record.tid).rsplit("|", 1)
        left, right = map(int, mark.split(":"))

        for pos in record.positions:
            base = record.query[pos - record.pos - record.qstart]
            pos += left  # Add offset of the aligned genomic region.

            if base not in "ACGT" or (chrom, pos) not in g:
                continue

            if getattr(g[chrom, pos], base) is not None:
                snps[chrom, pos][base] += record.mapq

    scores = itertools.chain(*(candidates.itervalues()
                               for candidates in snps.itervalues()))
    threshold = np.median(np.fromiter(scores, np.float64))

    # Filter the resulting mapping and output a set of the corresponding
    # genome identifiers.
    tsv = csv.writer(sys.stdout,
                     ["rsid", "chromosome", "position", "genotype"],
                     delimiter="\t")
    for (chrom, pos), candidates in snps.iteritems():
        if all(score <= threshold for score in candidates.itervalues()):
            continue

        # Order candidates by quality and pick a single allele in the
        # monoploid case, or two alleles in the diploid one.
        candidates = sorted(candidates.iteritems(),
                            key=lambda (_, quality): quality,
                            reverse=True)
        genotype = "".join(base for (base, _) in candidates[:diploid + 1])
        tsv.writerow([g[chrom, pos].id,
                      chrom, str(pos + 1), genotype])

    # Should we store frequency of an SNP in a particular sample in the
    # index?
    return snps


if __name__ == "__main__":
    globaloptions = [
        ("q", "quiet", False, "do not produce any output while running"),
        ("i", "index-root", None, "path to 'consequnce' index")
    ]
    opster.dispatch(globaloptions=globaloptions)
