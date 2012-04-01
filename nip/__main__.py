# -*- coding: utf-8 -*-

from __future__ import absolute_import, unicode_literals

import os.path
from collections import deque

import opster

import nip

#: Options, valid for each of `nip` commands.
globaloptions = [("q", "quiet", False, "suppress any possible output"),
                 ("c", "config", "config.json", "path to `nip` config file")]


@nip.setup
def index(config, read, *reads):
    def popread(sep="_"):
        """Returns the next 'valid' read in a sorted read list."""
        while reads:
            read = reads.popleft()
            fname, _ = os.path.splitext(read)
            if sep not in fname:
                continue  # Oops, bad file, captain?

            name, idx = fname.rsplit(sep, 1)
            if name in seen:
                continue  # We already have a pair for this read!

            if idx not in "12":
                continue  # Invalid pair index.

            return (name, idx), read
        else:
            raise IndexError

    reads = deque(sorted([read] + list(reads)))
    pairs, seen = [], set()
    while reads:
        try:
            (name1, idx1), read1 = popread()
            (name2, idx2), read2 = popread()

            if name1 == name2 and idx1 != idx2:
                pairs.append((read1, read2))
                seen.add(name1)
            else:
                reads.appendleft(read2)
        except IndexError:
            continue  # No moar reads?

    nip.index(config, pairs)


if __name__ == "__main__":
    opster.dispatch(globaloptions=globaloptions)
