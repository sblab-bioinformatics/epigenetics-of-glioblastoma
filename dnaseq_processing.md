<!-- MarkdownTOC -->

- [Software tools](#software-tools)
- [Alignment](#alignment)
- [Variation](#variation)
- [Tumour purity](#tumour-purity)
- [SNVs in ATRX gene](#snvs-in-atrx-gene)

<!-- /MarkdownTOC -->

Software tools
==============

Required software:

* [Isaac](http://bioinformatics.oxfordjournals.org/content/29/16/2041.long)
* [Strelka](http://bioinformatics.oxfordjournals.org/content/28/14/1811.long)
* [Canvas](http://bioinformatics.oxfordjournals.org/content/32/15/2375.long)
* [Manta](http://bioinformatics.oxfordjournals.org/content/32/8/1220.long)
* [SciClone](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003665)
* [Python](https://www.python.org/) version 2.7.10
* Human reference genome GRCh37 obtained from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html)


Alignment
=========

Quality control and alignment of reads from the HiSeq Illumina sequencer to the human reference genome GRCh37 were performed using `Isaac`. 


Variation
=========

<!--
This part from https://github.com/sblab-bioinformatics/projects/blob/b213a1517a54b9b7cd2726378faf2600897c0853/20150501_methylation_brain/20160202_snv_5hmC.md
-->

Somatic variant detection of the match tumour-margin sample obtained from the patient was performed using software developed by [Illumina](http://www.illumina.com/informatics/research.html). 

* Single nucleotide variations (SNVs) and small somatic indels were identified using `Strelka` from the aligned the sequencing reads.
* Large copy number variants (CNVs) and structural variants (SVs) were detected using `Canvas` and `Manta` respectively.

Changing the output of Strelka to `.bedGraph` format:

```python
import vcf # https://pyvcf.readthedocs.org/en/latest/index.html

vcf_reader = vcf.Reader(open('CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf', 'r')) # Output from Strelka

snv_tvm=[]

for record in vcf_reader:
    if record.is_snp:
        chr=record.CHROM
        start=record.start
        end=record.end
        base_m=record.REF
        base_t=str(record.ALT).replace("[","").replace("]","")
        qss=record.INFO["QSS"]
        tqss=record.INFO["TQSS"]
        snv_tvm.append((chr, start, end, base_m, base_t, qss, tqss))

len(snv_tvm) # 8169 SNVs in the tumour vs. margin comparison

bed_tvm=open("CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph", "w")

for snv in snv_tvm:
    bed_tvm.write("%s\n" % ("\t".join(["chr"+snv[0], str(snv[1]), str(snv[2]), snv[3], snv[4], str(snv[5]), str(snv[6])])))

bed_tvm.close()
```


Tumour purity
=============

First, the copy number variant analysis performed using `Canvas` provided us with an estimate of 71% purity of the tumour sample (PurityModelFit=0.0078 and InterModelDistance=0.0355). 

This software calculates coverage and tumour SNV allele frequencies at heterozygous germline positions along the genome. It then assigns copy numbers per genomic regions and infers genome-wide ploidy and purity by fitting the data to expected models for each copy number state given purity and diploid coverage level combinations. Purity is then derived from the best fitting purity/ploidy model. The lower PurityModelFit and InterModelDistance are, the better the fit of the data with the model. On truth data set, samples with PurityModelFit < 0.02 had >90% accurate CNV calls.

Second, purity was estimated from the somatic SNVs directly:

* High quality SNVs in obvious diploid regions (chromosomes 2,3,4,5,8,12,13,14,15,16,21) were selected.

* Subclones were investigated using `sciClone` to cluster the SNVs based on variant read frequency (VRF) and coverage depth (filtering anything below 15x). The software detected only one cluster, and removed variants detected as noise.

* Purity was estimate by multiplying by 2 the mean VRF of this cluster of high quality diploid SNVs (if the sample was pure, the mean would be 50%): 36.6*2 = 73.2%

Overall the tumour DNA shows a normal:tumour ratio of 30:70.


SNVs in ATRX gene
=================

<!--
This part from https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160826_genomic_profiles_examples/scripts/20160826_genomic_profiles_examples.md
-->

```R
library(data.table)

CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic <- fread("CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph")
setnames(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic, c("chr", "start", "end", "ref", "alt", "qss", "tqss"))

data.frame(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic[chr == "chrX" & start > 76500000 & end < 77400000])
# chrX  76799744  76799745   T   A  48    1 (intron)
# chrX  77030992  77030993   C   A  79    1 (intron)
```

