NGS Search
==========

Introduction
------------

Existing tools for SNP calling use Bayesian inference to pick the
*most likely* genotype for each particular site [2]; unfortunately, almost
all of these tools (GATK, samtools, SOAPsnp) require preprocessing of raw
sequencing data, to facilitate variant calling. This applies to both new
SNP discovery and finding already known SNPs.

In general, obtaining a complete and accurate variation record from
sequencing data can be done in three phases [1]:

* Phase 1. Raw reads are mapped to their correct origin in the reference
           genome (GATK authors mention 'elimination of molecular duplicates'
           followed by local realignment; this is only relevant to indel
           calling, right?)
* Phase 2. Resulting alignments are analyzed for statistically-significant variant
           sites, with alternate alleles present (SNP, indels, CNV).
* Phase 3. Raw variant calls from the previous step are then refined, using
           for example population-wide allele frequencies.

In this work we develop a method of calling SNPs directly from reads, bypassing
phases 2 and 3 of the above pipeline.

[1]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3083463/
[2]: http://www.nature.com/nrg/journal/v12/n6/full/nrg2986.html


Methods
-------

### Index organization

We maintain a separate *index* for each human chromosome. The index is
stored in a [persistent trie] [3], because:

* persistence allows for easy serialization and parallel independent
  updates for each of the new genomes;
* trie guarantees /O(bits of position)/ lookups and updates for each
  genomic position and reduced memory consumption.

Trie is indexed by variant sites' positions in the reference genome.
Each position is mapped to a 4-element tuple, representing possible
mutations for each of the genomic bases: `(A, C, G, T)`. One of the four
positions is **always** empty, because it conforms to the reference
base; the rest contain a possibly empty list of genome identifiers
with a corresponding variant. For example:

            A   C     G                 T
            |   |     |                 |
    141 -> ([], NULL, ["HTC10499_s_8"], [])
                ^
                |
           reference base

Here, we have a mapping for position `144`; reference base is `C` and a
single indexed genome with id `"HTC10499_s_8"` is present for variant
`C->G`.

[3]: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.5452

### Pipeline

1. Each read from the sequencing data is aligned to the *partial* reference
   genome (details on how this sequence is obtained are given later).
   No local realignment is necessary, since we only target SNPs and not
   indels or CNVs.
2. After the read is mapped, we know its exact position in the reference
   genome (including both chromosome name and base number). Thus, for each
   aligned read we lookup the position in the corresponding chromosomal
   index.

### Lookups

To lookup an aligned read (encoded in BAM or SAM) in the chromosomal index
we first unpack the CIGAR string, as described in
[SAM Format Specification] [4], ignoring everything but *Match-Mismatch*
blocks. This gives us a list of aligned positions in the corresponding
segment sequence. Then for each position we fetch a 4-element tuple from
the index and **if the base in the segment sequence matches any
non-reference bases** we cache it along with the looked up position and
the original read.

After all reads are processed in this way we have a cache of all the
possible SNPs, each of which has a different *quality*. SNP quality
is calculated from SNP frequency in the processed reads, aligment
quality (MAP!) for the reads cached with a particular SNP and nucleotide
coverage of the corresponding genomic position. So, all we have to
do is -- pick for each genomic position an allele with the highest
*quality*.

#### A note on diploidy

The lookup procedure described above can be easily extended to handle
diploid genomes. We just need to pick *two* alleles with the highest
*quality* for each genomic position, instead of a single allele.

#### A note on *partial* reference

Because we use a *partial* reference sequence, aligned positions in the
SAM or BAM file should be converted, to the corresponding genomic
positions. This can be easily done, by embedding the corresponding
genomic region into sequence identifier, for example:

    >gi|49175990|ref|NC_000913.2||222716:224764
                                  ^      ^
                                  |      |
                                start   end

[4]: http://samtools.sourceforge.net/SAM1.pdf

### Updating the index

#### Adding new genomes

Updating the index with new SNPs can be done in the background, following
the full pipeline, described above. Because the index structure is
*persistent*, the update can be done in parallel, independently for each
new genome to be indexed.

#### Rebuilding *partial* reference

To speedup the alignment step we only store meaningful regions of the
reference genome; that is -- regions with at least one indexed variant
per doubled insert size (strictly speaking, the insert size varies from
platform to platform and from run to run, so we use double maximum
insert size, found in the latest 1000genomes [release] [5] -- **7000**
base pairs).

[5]: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/sequence.index

Results
-------

* Compare to [Crossbow] [6] and [SNiPlay] [7].

[6]: http://www.ncbi.nlm.nih.gov/pubmed/19930550
[7]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3102043/
