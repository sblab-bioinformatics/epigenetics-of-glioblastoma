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
* [SciClone](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003665)
* Human reference genome GRCh37 obtained from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html)


Alignment
=========

Quality control and alignment of reads from the HiSeq Illumina sequencer to the human reference genome GRCh37 were performed using `Isaac`. 


Variation
=========

Somatic variant detection of the match tumour-margin sample obtained from the patient was performed using software developed by [Illumina](http://www.illumina.com/informatics/research.html). 

* Single nucleotide variations (SNVs) and small somatic indels were identified using `Strelka` from the aligned the sequencing reads.
* Large copy number variants (CNVs) and structural variants (SVs) were detected using `Canvas` and `Manta` respectively.


Tumour purity
=============

First, the copy number variant analysis performed using `Canvas` provided us with an estimate of 71% purity of the tumour sample (PurityModelFit=0.0078 and InterModelDistance=0.0355). 

This software calculates coverage and tumour SNV allele frequencies at heterozygous germline positions along the genome. It then assigns copy numbers per genomic regions and infers genome-wide ploidy and purity by fitting the data to expected models for each copy number state given purity and diploid coverage level combinations. Purity is then derived from the best fitting purity/ploidy model. The lower PurityModelFit and InterModelDistance are, the better the fit of the data with the model. On truth data set, samples with PurityModelFit < 0.02 had >90% accurate CNV calls.

Second, purity was estimated from the somatic SNVs directly:

* High quality SNVs in obvious diploid regions (chromosomes 2,3,4,5,8,12,13,14,15,16,21) were selected.

* Subclones were investigated using `sciClone` to cluster the SNVs based on variant read frequency (VRF) and coverage depth (filtering anything below 15x). The software detected only one cluster, and removed variants detected as noise.

* Purity was estimate by multiplying by 2 the mean VRF of this cluster of high quality diploid SNVs (if the sample was pure, the mean would be 50%): 36.6*2 = 73.2%

Overall the tumour DNA shows a normal:tumour ratio of 30:70.
