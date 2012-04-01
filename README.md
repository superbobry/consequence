ngs-search
==========

Draft variance calling pipeline
-------------------------------

Useful BioStar Threads:

* [What is the best pipeline for human whole exome sequencing?](http://biostar.stackexchange.com/questions/1269/what-is-the-best-pipeline-for-human-whole-exome-sequencing)
* [workflow or tutorial for SNP calling?](http://biostar.stackexchange.com/questions/8260/workflow-or-tutorial-for-snp-calling)

Available options:

* [pipette](https://github.com/metalhelix/pipette) -- Ruby framework for
  variance calling pipelines

Current vision:

1. Align reads to the reference
   [Grch37](http://www.1000genomes.org/category/frequently-asked-questions/grch37)
   using Bowtie; this can be done in parallel for each of the available SRR's.
2. Convert the resulting SAM file to a indexed & sorted BAM:

        $ samtools view -bS <input.sam> > output.bam
        $ samtools sort <input.bam> <input.sorted>
        $ samtools index <input.sorted.bam>

3. (optional, since we don't target indels at the moment) Perform local
   re-alignment using `LocalRealigner` from
   [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels)
   and [call](http://www.broadinstitute.org/gsa/wiki/index.php/Indel_Genotyper_V2.0)
   indels.
4. Call SNPs with either
   [UnifiedGenotyper](http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper)
   or `samtools mpileup`
