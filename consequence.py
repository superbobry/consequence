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
import operator
import os.path
import sys
from collections import defaultdict, namedtuple, MutableMapping

import vcf
import pysam
from Bio import SeqIO


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
            path = os.path.join(self.base_path, "{0}.index".format(chr))
            cPickle.dump(self.cache[chr], open(path, "wb"))

    def load(self, *chrs):
        if not chrs:
            fnames = glob.glob(os.path.join(self.base_path, "*.index"))
            chrs = [os.path.basename(name)
                    for name, _ext in map(os.path.splitext, fnames)]

        for chr in chrs:
            if chr in self.cache:
                continue  # Oops, already got that sequence loaded.

            path = os.path.join(self.base_path, "{0}.index".format(chr))
            if os.path.isfile(path):
                # Note(Sergei): this is probably the dumbest thing possible,
                # but hey -- it's still better than nothing ;)
                self.cache[chr] = cPickle.load(open(path, "rb"))
            else:
                # Yeah-yeah, it's not a 'persistent trie', because there
                # isn't one avaiable for Python yet.
                self.cache[chr] = {}

        return [self.cache[chrom] for chrom in chroms]

    def __repr__(self):
        return "<Genome: {0!r}>".format(self.cache)

    def __len__(self):
        # Obviously, this does *NOT* include non-cached indices.
        return sum(map(len, self.cache))

    def __iter__(self):
        return ((chr, pos) for pos in self.cache[chr] for chr in self.cache)

    def __contains__(self, (chr, pos)):
        if chr not in self.cache:
            self.load(chr)

        return pos in self.cache[chr]

    def __getitem__(self, (chr, pos)):
        if chr not in self.cache:
            self.load(chr)

        return self.cache[chr][pos]

    def __setitem__(self, (chr, pos), chunk):
        if chr not in self.cache:
            self.load(chr)

        self.cache[chr][pos] = chunk

    def __delitem__(self, (chr , pos)):
        del self.cache[chr][pos]


def update_index(seq_path, seq_id):
    g = Genome()

    for record in vcf.VCFReader(open(seq_path, "rb")):
        if not record.is_snp:
            continue

        # XXX here comes the punchline, VCF uses 1-based indexing,
        # so does SAM / BAM, but 'pysam' converts *all* indexes to
        # 0-based!
        chrom, pos = record.CHROM, record.POS - 1
        g.setdefault((chrom, pos), Chunk(**{record.REF: False}))

        for alt in filter(operator.truth, record.ALT):
            known = getattr(g[chrom, pos], alt) or set()
            if seq_id in known:
                continue

            known.add(seq_id)
            g[chrom, pos] = g[chrom, pos]._replace(**{alt: known})

    g.dump()


def build_partial_reference(ref_path, insert_size=1024):
    if not os.path.isfile(ref_path):
        raise RuntimeError("Invalid reference sequence: {0!r}"
                           .format(ref_path))

    g = Genome()
    cropped = []
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


def naive_lookup(seq_path, is_diploid=False):
    g = Genome()
    f = pysam.Samfile(seq_path, "rb")

    cov  = defaultdict(lambda: defaultdict(int))
    snps = defaultdict(lambda: defaultdict(int))
    for record in f:
        if record.is_duplicate or record.is_unmapped or record.tid < 0:
            continue

        chrom = f.getrname(record.tid)

        for pos in record.positions:
            base = record.query[pos - record.pos - record.qstart]

            if base not in "ACGT" or (chrom, pos) not in g:
                continue

            if getattr(g[chrom, pos], base) is not None:
                cov[pos][base] += 1
                snps[chrom, pos][base] += record.mapq

    # Calculate scores for each possible SNPs.
    threshold = 0.
    for (chrom, pos), candidates in snps.iteritems():
        for base in candidates:
            # Resulting qual. score is proportional to the base coverage.
            candidates[base] /= float(cov[pos][base])
            threshold = max(candidates[base], threshold)

    threshold /= 2.

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
        genotype = "".join(base for (base, _) in candidates[:is_diploid + 1])
        tsv.writerow(["*", chrom, str(pos), genotype])

    # Should we store frequency of an SNP in a particular sample in the
    # index?
    return snps


if __name__ == "__main__":
    # build_partial_reference(ref_seq)
    # update_index(seq_path, seq_id)
    # naive_lookup(seq_path)
    pass
