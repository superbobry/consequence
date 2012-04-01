# -*- coding: utf-8 -*-
"""
    nip
    ~~~

    See 'README.md' for details.
"""

from __future__ import unicode_literals

__all__ = ["index"]

import collections
import json
import os.path
import tempfile
from functools import wraps

import logbook
import opster
import pysam
from iterpipes import cmd, check_call
from logbook import StderrHandler, NullHandler


def index(config, reads):
    """Index a given list of paired ``reads`` relative to a given
    ``reference`` genome.

    :param Config config: current config object.
    :param list reads: a list of paired reads, where each element is a
                       tuple of two elements, for example:
                       ``("foo_1.fastq", "foo_2.fastq")``.
    """
    # For each read:
    # * align to the reference
    # * sort and index the resulting BAM file
    sorted_alignments = []

    for alignment in [align_reads(config, left, right) for left, right in reads]:
        fname = tempfile.mktemp(suffix=".bam", dir=config["tmp_dir"])
        pysam.sort(alignment, fname)
        os.remove(alignment)
        pysam.index(fname)
        sorted_alignments.append(fname)

    # Call SNPs in the sorted BAM files.
    variances = call_snps(config, sorted_alignments)

    # TODO(Sergei): filter low quality SNPs?
    # TODO(Sergei): merge resulting VCF file with 'nip.index'?


# Internal.

class Config(collections.Mapping):
    """A minimal wrapper around JSON config file."""

    def __init__(self, defaults=()):
        self.data = dict(defaults)

    def __repr__(self):
        return "<Config: {0!r}>".format(self.data)

    def __str__(self):
        return "<Config: {0}>".format(self.data)

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, key):
        return self.data[key]

    @classmethod
    def from_file(cls, path):
        if os.path.exists(path):
            with open(path) as f:
                return cls(json.load(f))
        else:
            raise RuntimeError("Config file {0} doesn't exist!'")


def setup(f):
    """Setup a given `nip` command."""
    @wraps(f)
    def inner(*args, **kwargs):
        config = Config.from_file(kwargs.pop("config"))

        if not kwargs.pop("quiet"):
            handler = StderrHandler(level="INFO")
        else:
            handler = NullHandler()

        with handler.applicationbound():
            f(config, *args)

    return opster.command()(inner)


def align_reads(config, left, right):
    """Run `bowtie` on a given read pair and returns a path to the
    resulting BAM file.
    """
    logbook.info("Aligning [{0}, {1}]", left, right)
    fname = tempfile.mktemp(suffix=".bam", dir=config["tmp_dir"])
    check_call(cmd("bowtie "
                   "{0[alignment_options]} "
                   "{0[reference_index]} "
                   "-1 {1} -2 {2} --sam | samtools view -bS - > {3}"
                   .format(config, left, right, fname)))
    return fname


def call_snps(config, sorted_alignments):
    """Run `UnifiedGenotyper` on a given list of sorted BAM files."""
    logbook.info("Calling SNPs on {0}", ", ".join(sorted_alignments))
    fname = tempfile.mktemp(suffix=".vcf", dir=config["tmp_dir"])
    check_call(cmd("GenomeAnalysisTK.jar "
                   "-R {0[reference_fasta]} "
                   "-T UnifiedGenotyper "
                   "--dbsnp {0[dbsnp]} ",
                   "-o {1}".format(config, fname) +
                   " ".join("-I {0}".format(a) for a in sorted_alignments)))
    return fname
