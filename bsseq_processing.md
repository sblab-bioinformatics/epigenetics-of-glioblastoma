<!-- MarkdownTOC -->

- [Software tools and requirements](#software-tools-and-requirements)
- [Read trimming and alignment](#read-trimming-and-alignment)
- [Merging bam files and marking duplicates](#merging-bam-files-and-marking-duplicates)
- [Methylation at CpG sites](#methylation-at-cpg-sites)
- [TODO](#todo)

<!-- /MarkdownTOC -->

Software tools and requirements
===============================

All data processing was implemented under Linux environment with standard GNU coreutils tools. 
Additional required software:

* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html) version 1.9
* [bwa-meth](https://github.com/brentp/bwa-meth)
* [bwa](https://github.com/lh3/bwa) version 0.7
* [samtools](http://www.htslib.org/) version 1.1
* [picard/MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html) version 1.127
* [bam2methylation.py](https://github.com/dariober/bioinformatics-cafe/tree/master/bam2methylation)
* [BamUtil: clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)
* [bedtools](http://bedtools.readthedocs.org/en/latest/) version 2.25
* [tabix](http://www.htslib.org/doc/tabix.html) version 0.2.5
* [R-3.2.3](https://cran.r-project.org/) 
* Human reference genome version hg19 obtained from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

Read trimming and alignment
===========================

<!--
This part from 
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151208_BSnewdata/20151208_BSnewdata.md
-->

Raw reads from the HiSeq Illumina sequencer were trimmed to remove adapters ligated to the 3'end 
with `cutadapt` and aligned to the human reference genome version hg19 with `bwameth`. After alignment, the overlap between read pairs 
was soft-clipped using `bam clipOverlap` in order to avoid double-counting fragments sequenced twice by the same read pair. 

The reference genome was first indexed for bwameth with:

```
bwameth.py index genome.fa
```

Then, reads were trimmed, aligned and the overlap clipped:

```
cutadapt -m 10 -O 1 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o fastq_trimmed/$out1 -p fastq_trimmed/$out2 $fq1 $fq2 && 
bwameth.py -t 4 --reference genome.fa --prefix bwameth/20151208_BSnewdata/${out}.unclipped fastq_trimmed/$out1 fastq_trimmed/$out2 && 
bam clipOverlap --in bwameth/20151208_BSnewdata/${out}.unclipped.bam --out bwameth/20151208_BSnewdata/${out}.bam --stats --storeOrig XC && 
```

`$fq1` and `$fq2` are variables for fastq file mate 1 and mate 2, respectively. Similarly, `$out1` and `$out2` are variables for the 
name of the trimmed fastq files. `$out` is the basename for the aligned output files.

Merging bam files and marking duplicates
========================================

<!--
This part from 
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151208_BSnewdata/20151208_BSnewdata.md
-->


Bam files generated from the same library but sequenced on different lanes were merged and positional duplicates were marked:

```
samtools merge -f -@8 ${id}.${runBatch}.bam $bams &&
samtools index ${id}.${runBatch}.bam &&
java -Xmx3G -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT TMP_DIR=./ I=${id}.${runBatch}.bam O=/dev/null M=${id}.${runBatch}.markdup.txt
```

`${id}.${runBatch}` is the basename for the merged and marked bam files. `$bams` is the list of bam files to be merged.


Methylation at CpG sites
========================

<!--
This part from 
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151208_BSnewdata/20151208_BSnewdata.md
-->


Methylation at each CpG site was extracted from bam files using `bam2methylation.py`. 
Reads were included in the methylation calling if their mapping quality was 10 or above 
(mapping quality 10 corresponds to a probability of being wrongly mapped of 10%). Alignments were excluded if
the aligner marked them as unmapped, or not primary, or failed quality check, or supplementary (flag 2820).
In addition read bases with quality below 15 were excluded.

```
## Process each chromsome separately
chroms=`cut -f 1 genome.fa.fai`

mkdir chroms
for bam in *.bam
do
    for chr in $chroms
    do
        bdg=${bam%.bam}.$chr.bedGraph
        bam2methylation.py -l hg19.allCpG.bed.gz -i $bam -r genome.fa -A -s ' -q 10 -F2820' -mm -mq 15 -R $chr \
        | pigz -f > chroms/${bdg}.gz" 
    done
done
```

`hg19.allCpG.bed.gz` is a bed file containing the position of all the CpG in the reference genome.


Concatenate files within libraries and index:

```
zcat ${id}.chr*.bedGraph.gz | bgzip > ${id}.cpg.bedGraph.gz &&
tabix -p bed ${id}.cpg.bedGraph.gz"
```

`${id}` is the library basename


TODO
====

Code for all figures and tables, depending on what we put in the manuscript
