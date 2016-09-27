<!-- MarkdownTOC -->

- [Software tools](#software-tools)
- [Alignment](#alignment)
- [Variation](#variation)
- [Tumour purity](#tumour-purity)

<!-- /MarkdownTOC -->

Software tools
==============

Required software:

* [Isaac](http://bioinformatics.oxfordjournals.org/content/29/16/2041.long)
* [Strelka](http://bioinformatics.oxfordjournals.org/content/28/14/1811.long)
* [Canvas](http://bioinformatics.oxfordjournals.org/content/32/15/2375.long)
* [Manta](http://bioinformatics.oxfordjournals.org/content/32/8/1220.long)
* Human reference genome GRCh37 obtained from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html)


Alignment
=========

Quality control and alignment of reads from the HiSeq Illumina sequencer to the human reference genome GRCh37 were performed using `Isaac`. 


Variation
=========

Somatic variant detection of the match tumour-margin sample obtained from the patient was performed using software developed by [Illumina](http://www.illumina.com/informatics/research.html). Single nucleotide variations (SNVs) and small somatic indels were identified using `Strelka` from the aligned the sequencing reads. Large copy number variants (CNVs) and structural variants (SVs) were detected using `Canvas` and `Manta` respectively.


Tumour purity
=============

Check email from Jenn.

