# -*- coding: utf-8 -*-
"""
    nip
    ~~~

    See 'README.md' for details.
"""

from __future__ import unicode_literals

__all__ = ["index"]


def index(reference, reads):
    """Index a given list of paired ``reads`` relative to a given
    ``reference`` genome.

    :param str reference: absolute or relative path to the reference
                          FASTA file.
    :param list reads: a list of paired reads, where each element is a
                       tuple of two elements, for example:
                       ``("foo_1.fastq", "foo_2.fastq")``.
    """
