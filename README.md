NGS Search
==========

Note: '*' выделены места, которые еще предстоит дописать.

Introduction
------------

* Importance of SNP calling, applications;
* Major problems for accurate SNP calling; current trend -- a lot of cheap
  but inaccurate reads, which leads to unavoidable errors on all phases of
  variant detection;

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

We maintain a separate *index* for each human chromosome. An index is
basically a [persistent trie] [3], indexed by variant sites' positions
in the reference genome. Each position is mapped to a 4-element tuple,
representing possible mutations for each of the genomic bases: `(A, C, G, T)`.
One of the four positions is **always** empty, because it conforms
to the reference base; the rest contain a possibly empty list of
genome identifiers with a corresponding variant. For example:

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

* Describe how low coverage variants are filtered out.

### Updating the index

#### Adding new genomes

Updating the index with new SNPs can be done in the background, following
the full pipeline, described above. Because the index structure is
*persistent*, the update can be done in parallel, independently for each
new genome to be indexed.

#### Rebuilding *partial* reference

To speedup the alignment step we only store meaningful regions of the
reference genome; that is -- regions with at least one indexed variant
per doubled read length (150 base pairs for Illumina reads?). An obvious
downside of this approach is that resulting alignments will have positions
relative to the *partial* reference, instead of the original one. So an
**extra index**, mapping intervals from *partial* reference to the original
reference, is required.

### Accuracy

* Is it possible to estimate accuracy of our approach? on first thought
  it follows directly from reads quality and alignment accuracy (see also
  note above the quality of NGS reads);

Results
-------

* Compare to Crossbow [4] and SNiPlay [5].

[4]: http://www.ncbi.nlm.nih.gov/pubmed/19930550
[5]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3102043/
