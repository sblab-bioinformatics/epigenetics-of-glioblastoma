<!-- MarkdownTOC -->

- [Software tools and requirements](#software-tools-and-requirements)
- [Transcript expression from RNA-Seq](#transcript-expression-from-rna-seq)
- [TODO](#todo)

<!-- /MarkdownTOC -->

Software tools and requirements
===============================

All data processing was implemented under Linux environment with standard GNU coreutils tools. 
Additional required software:

* [kallisto](https://pachterlab.github.io/kallisto/)
* [R-3.2.3](https://cran.r-project.org/) 
* Transcripts in fasta format for [hg19/grch37](ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz)

Transcript expression from RNA-Seq
==================================

<!-- 
For this part see 
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160303_rnaseq/20160303_rnaseq.md
 -->

Transcript expressions were assessed by aligning RNA-Seq against the reference transcripts from Ensembl using `Kallisto`. 
First reference sequences were indexed:

```
kallisto index -i Homo_sapiens.GRCh37.cdna.all.idx Homo_sapiens.GRCh37.cdna.all.fa.gz
```

Transcript quantitation was obtained for each library as follows:

```
kallisto quant --bias -i Homo_sapiens.GRCh37.cdna.all.idx -o ${x}.tx <(zcat $fq1) <(zcat $fq2)
```

Where `$fq1` and `$fq2` are fastq file for the two mates of each pair and `${x}` is the base-name for the output file.

TODO
====

Code for all figures and tables, depending on what we put in the manuscript