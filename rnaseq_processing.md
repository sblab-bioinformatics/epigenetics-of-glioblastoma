<!-- MarkdownTOC -->

- [Software tools and requirements](#software-tools-and-requirements)
- [Transcript expression from RNA-Seq](#transcript-expression-from-rna-seq)
- [Transcript levels ATRX gene](#transcript-levels-atrx-gene)
- [TODO](#todo)

<!-- /MarkdownTOC -->

Software tools and requirements
===============================

All data processing was implemented under Linux environment with standard GNU coreutils tools. 
Additional required software:

* [kallisto](https://pachterlab.github.io/kallisto/)
* [R-3.2.3](https://cran.r-project.org/) 
* Transcripts in fasta format for hg19/grch37(ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz)
* Annotation [mart_export.GRCh37.txt](http://grch37.ensembl.org/biomart/martview/) to link ensembl transcripts to gene name with the appropropriate attributes. Header line edited to replace whitespace with underscore and parenthesis remove.

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

To generate a table of transcripts quantification:

```R
library(data.table)
library(ggplot2)

tx<- fread('tableCat.py -H -i ear04*.tx_quant.tsv -r "\\..*" -id "library_id"') # ear04*.tx_quant.tsv are obtained from the kallisto alignment after renaming
ann<- fread('zcat mart_export.GRCh37.txt.gz')
tx[, Ensembl_Transcript_ID := sub('\\..*', '', target_id)]
tx<- merge(tx, ann[], by= 'Ensembl_Transcript_ID')

## Keep only regular chroms and rename to UCSC format. 
## Note this will make the sum(TPM) < 1e6
tx<- tx[Chromosome_Name %in% c(1:22, 'MT', 'X', 'Y')]
tx[, chrom := ifelse(Chromosome_Name %in% c(1:22, 'X', 'Y'), paste0('chr', Chromosome_Name), Chromosome_Name)]
tx[ , chrom := ifelse(Chromosome_Name == 'MT', 'chrM', chrom)]
stopifnot(unique(sort(tx$chrom)) == c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrM", "chrX", "chrY"))
tx[, Strand := ifelse(Strand == 1, '+', '-')]

## Note Change to bed coords!!
tx[, Transcript_Start_bp := Transcript_Start_bp-1]
tx[, Gene_Start_bp := Gene_Start_bp-1] 

## Table logFC
## ===========
txLfc<- dcast.data.table(data= tx[library_id %in% c("ear047_F3", "ear049_M3")], Ensembl_Transcript_ID ~ library_id, value.var= 'tpm')
## Note using unique(..., by= NULL)!!!
txLfc<- merge(txLfc, 
            unique(tx[, list(Ensembl_Transcript_ID, chrom, Transcript_Start_bp, Transcript_End_bp, Associated_Gene_Name, Strand)], by= NULL), 
        by= 'Ensembl_Transcript_ID')

txLfc[, tpmLog2FC := log2((ear047_F3 + 0.01) / (ear049_M3 + 0.01))]
txLfc$tpmLog2Avg<- log2(rowMeans(txLfc[, list(ear047_F3, ear049_M3)]) + 0.01)
txLfc[, promStart := ifelse(Strand == '+', Transcript_Start_bp-1000, Transcript_End_bp)]
txLfc[, promEnd := promStart + 1000]
write.table(file= 'tx_quant_lfc.bed',
    txLfc[, list(chrom, promStart, promEnd, Ensembl_Transcript_ID, Associated_Gene_Name, Strand, Transcript_Start_bp, Transcript_End_bp, tpmLog2FC, tpmLog2Avg, ear047_F3, ear049_M3)][order(chrom, Transcript_Start_bp, Transcript_End_bp)],
    row.names= FALSE, quote= FALSE, sep= '\t')
```


Transcript levels ATRX gene
===========================

<!-- 
From https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160826_genomic_profiles_examples/scripts/20160826_genomic_profiles_examples.md
 -->

```R
library(data.table)
library(ggplot2)

tx<- fread('tx_quant_lfc.bed')

gxPlus<- tx[Strand == '+', list(
    chrom= min(chrom),
    promStart= min(promStart),
    Strand= min(Strand),
    ear047_F3= sum(ear047_F3),
    ear049_M3= sum(ear049_M3)), by= list(Associated_Gene_Name)]

gxMinus<- tx[Strand == '-', list(
    chrom= min(chrom),
    promStart= max(promStart),
    Strand= min(Strand),
    ear047_F3= sum(ear047_F3),
    ear049_M3= sum(ear049_M3)), by= list(Associated_Gene_Name)]

gx<- rbindlist(list(gxPlus, gxMinus))
gx[, promEnd := promStart + 1000]
gx[, tpmLog2FC := log2((ear047_F3 + 0.01) / (ear049_M3 + 0.01))]
gx$tpmLog2Avg<- rowMeans(gx[, list(log2(ear047_F3 + 0.01), log2(ear049_M3 + 0.01))])
gx<- gx[, list(chrom, promStart, promEnd, Associated_Gene_Name, Strand, ear047_F3, ear049_M3, tpmLog2FC, tpmLog2Avg)][order(chrom, promStart, promEnd)]


xdata2 <- gx[Associated_Gene_Name == "ATRX"]

tab2 <- xdata2[, list(Associated_Gene_Name = Associated_Gene_Name, ear047_F3.log = log2(ear047_F3 + 0.01), ear049_M3.log = log2(ear049_M3 + 0.01))]

tab.summary2 <- reshape2::melt(tab2, id.vars = "Associated_Gene_Name", variable.name = "sample", value.name = "expression")
tab.summary2[, sample2 := factor(c("Tumour", "Margin"))]
tab.summary2[, sample := NULL]

gg <- ggplot(tab.summary2, aes(x = sample2, y = expression)) +
geom_bar(stat = "identity") +
xlab("") +
ylab(expression("log"[2]*"(TPM+0.01)")) +
theme_classic() +
theme(legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_text(size=10), axis.title.y = element_text(size=10)) +
coord_cartesian(ylim = c(3, 6)) +
scale_y_continuous(breaks = seq(0, 10, 1))
ggsave("figure/ATRX.expression.ylim3_6.png", width = 4.25/2.54, height = 4.75/2.54)
```

<img src="figures/ATRX.expression.ylim3_6.png" width="300">


TODO
====

Code for all figures and tables, depending on what we put in the manuscript
