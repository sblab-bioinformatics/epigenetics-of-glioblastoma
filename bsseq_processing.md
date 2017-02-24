<!-- MarkdownTOC -->

- [Software tools and requirements](#software-tools-and-requirements)
- [Read trimming and alignment](#read-trimming-and-alignment)
- [Merging bam files and marking duplicates](#merging-bam-files-and-marking-duplicates)
- [Methylation at CpG sites](#methylation-at-cpg-sites)
- [5mC and 5hmC levels ATRX gene](#5mc-and-5hmc-levels-atrx-gene)
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
* [COSMIC](http://cancer.sanger.ac.uk/cosmic) version 76


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
tabix -p bed ${id}.cpg.bedGraph.gz
```

`${id}` is the library basename.


5mC and 5hmC levels ATRX gene
=============================

<!--
This part modified from https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160826_genomic_profiles_examples/scripts/20160826_genomic_profiles_examples.md and https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160226_tsg_ocg.md
-->

Plot ATRX gene:

```R
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)

## Define genomic coordinates track
gtrack <- GenomeAxisTrack()
gtrack


## Define chromosome X ideogram track
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
itrack


## Load gene models using UcscTrack. Coordinates from https://www.ncbi.nlm.nih.gov/gene?cmd=retrieve&dopt=default&rn=1&list_uids=546
from <- 76500000
to <- 77400000
refGenes <- UcscTrack(track = "refGene", trackType = "GeneRegionTrack", genome = "hg19", chromosome = "chrX", name = "Genes", from = from, to = to, rstarts = "exonStarts", rends = "exonEnds", gene = "name2", symbol = "name2", transcript = "name", strand = "strand")


## Plot
displayPars(itrack) <- list(background.panel = "white", showBandId = FALSE, fontcolor = "black", fontsize = 20)
displayPars(gtrack) <- list(background.panel = "white", add53=TRUE, add35=TRUE, labelPos = "above", col = "black", fontcolor = "black", fontsize = 25)
displayPars(refGenes) <- list(transcriptAnnotation = "symbol", collapseTranscripts = "longest", shape = "arrow", fontfamily = "sans", fontsize.group = 30, background.panel = "white", col = NULL, background.title = "white", fontcolor.group = "black", col.line = "black")

png("figures/ATRX.gene.png", width = 40/2.54, height = 5/2.54)
plotTracks(list(itrack, gtrack, refGenes), from = 76740000, to = 77052000)
dev.off()
```

<img src="figures/ATRX.gene.png" width="600">

Base resolution 5mC and 5hmC levels at CpGi in 5'-UTR:

```bash
for bed in bedGraphs
do
    bname=`cut -d"." -f1 $bed`
    echo -e "chrX\t77041003\t77041725" | bedtools intersect -a $bed -b - -wa > $bname.cpg.ATRX.CGI.bedGraph
done
```

`$bedGraphs` is a list of files containing the position of all the CpG in the reference genome and the methylation counts (ear042_M8BS.cpg.bedGraph, ear043_M8oxBS.cpg.bedGraph, ear044_T3BS.cpg.bedGraph and ear045_T3oxBS.cpg.bedGraph).

```R
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(Biostrings)

ear042_M8BS <- fread("ear042_M8BS.cpg.ATRX.CGI.bedGraph")
setnames(ear042_M8BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear043_M8oxBS <- fread("ear043_M8oxBS.cpg.ATRX.CGI.bedGraph")
setnames(ear043_M8oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear044_T3BS <- fread("ear044_T3BS.cpg.ATRX.CGI.bedGraph")
setnames(ear044_T3BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear045_T3oxBS <- fread("ear045_T3oxBS.cpg.ATRX.CGI.bedGraph")
setnames(ear045_T3oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

setkey(ear042_M8BS, chr, start, end, strand)
setkey(ear043_M8oxBS, chr, start, end, strand)
setkey(ear044_T3BS, chr, start, end, strand)
setkey(ear045_T3oxBS, chr, start, end, strand)

M8 <- merge(ear042_M8BS, ear043_M8oxBS, suffixes = c(".BS", ".oxBS"))
T3 <- merge(ear044_T3BS, ear045_T3oxBS, suffixes = c(".BS", ".oxBS"))

setkey(M8, chr, start, end, strand)
setkey(T3, chr, start, end, strand)
M8.T3 <- merge(M8, T3, suffixes = c(".M8", ".T3"))

strack <- SequenceTrack(Hsapiens, chromosome = "chrX")
as.character(subseq(strack, 77041004, 77041725)) # 77041003 + 1

nucleotides <- unlist(strsplit(as.character(subseq(strack, 77041004, 77041725)), split = ""))
chroms <- rep("chrX", length(nucleotides))
starts <- seq(77041003, 77041724)
ends <- seq(77041004, 77041725)

cgi <- data.table(data.frame(chr = chroms, start = starts, end = ends, base = nucleotides))

tab5 <- M8.T3[, list(pct_5mC.M8=100*sum(cnt_met.oxBS.M8)/sum(cnt_tot.oxBS.M8), pct_5hmC.M8=100*(sum(cnt_met.BS.M8)/sum(cnt_tot.BS.M8) - sum(cnt_met.oxBS.M8)/sum(cnt_tot.oxBS.M8)), pct_C_other.M8=100*(1-sum(cnt_met.BS.M8)/sum(cnt_tot.BS.M8)), pct_5mC.T3=100*sum(cnt_met.oxBS.T3)/sum(cnt_tot.oxBS.T3), pct_5hmC.T3=100*(sum(cnt_met.BS.T3)/sum(cnt_tot.BS.T3) - sum(cnt_met.oxBS.T3)/sum(cnt_tot.oxBS.T3)), pct_C_other.T3=100*(1-sum(cnt_met.BS.T3)/sum(cnt_tot.BS.T3))), by=list(chr, start, end)]

setkey(cgi, chr, start, end)
setkey(tab5, chr, start, end)

cgi_tab5 <- data.frame(merge(cgi, tab5, all = TRUE))

cgi_tab5[is.na(cgi_tab5)] <- 0 # change all NAs to 0s

cgi_tab5 <- data.table(cgi_tab5) # change back to data.table

cgi_tabsummary5 <- reshape2::melt(cgi_tab5, id.vars = c("chr", "start", "end", "base"), variable.name = "sample", value.name = "pct_met")
cgi_tabsummary5[, sample2 := factor(c(rep("M", 722*3), rep("T", 722*3)))]
cgi_tabsummary5[, label := factor(rep(c(rep("5mC", 722), rep("5hmC", 722), rep("C", 722)),2), levels = c("C", "5hmC", "5mC"))]
cgi_tabsummary5[, sample:=NULL]
cgi_tabsummary5[, index := rep(1:722, 6)]

vals <- paste(nucleotides, as.character(complement(DNAStringSet(nucleotides))), sep = "\n")
keys <- as.character(1:722)
names(vals) <- keys
labels <- vals

gg <- ggplot(cgi_tabsummary5[index > 525 & index < 558], aes(x = sample2, y = pct_met, fill = label)) +
geom_bar(stat = "identity") +
xlab("") +
ylab("% Modification") +
theme_classic() +
theme(legend.title = element_blank(), strip.background = element_blank(), strip.text = element_text(size=8), panel.margin = unit(0.05, "lines"), axis.text.x = element_text(size=6)) +
scale_fill_manual(values=c("lightgray", "#009E73", "#D55E00")) +
facet_grid(. ~ index, switch = "x", labeller = labeller(index = labels))
ggsave("figures/ATRX.cpg.CGI.BS.oxBS.chrX.77041528.77041560.png", width = 15/2.54, height = 5/2.54, limitsize = FALSE)
```

<img src="figures/ATRX.cpg.CGI.BS.oxBS.chrX.77041528.77041560.png" width="600">

Mean 5mC and 5hmC levels in promoter:

```bash
for bed in bedGraphs
do
    bname=`cut -d"." -f1 $bed`
    bedtools intersect -a $bed -b $bed_cancer_genes_promoters_1kb -wa -wb | cut -f 1,2,3,4,5,6,10,11,12 > $bname.cpg_genes.promoters.all.cancer.glioblastoma.1000.sorted.bed
done
```

where `$bed_cancer_genes_promoters_1kb` contains the coordinates of all promoters of cancer genes as defined in COSMIC version 76.

```R
library(data.table)
library(ggplot2)

ear042_M8BS <- fread("ear042_M8BS.cpg_genes.promoters.all.cancer.glioblastoma.1000.sorted.bed")
setnames(ear042_M8BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "gene_name", "cancer", "glioblastoma"))

ear043_M8oxBS <- fread("ear043_M8oxBS.cpg_genes.promoters.all.cancer.glioblastoma.1000.sorted.bed")
setnames(ear043_M8oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "gene_name", "cancer", "glioblastoma"))

ear044_T3BS <- fread("ear044_T3BS.cpg_genes.promoters.all.cancer.glioblastoma.1000.sorted.bed")
setnames(ear044_T3BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "gene_name", "cancer", "glioblastoma"))

ear045_T3oxBS <- fread("ear045_T3oxBS.cpg_genes.promoters.all.cancer.glioblastoma.1000.sorted.bed")
setnames(ear045_T3oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "gene_name", "cancer", "glioblastoma"))

setkey(ear042_M8BS, chr, start, end, strand, gene_name, cancer, glioblastoma)
setkey(ear043_M8oxBS, chr, start, end, strand, gene_name, cancer, glioblastoma)
setkey(ear044_T3BS, chr, start, end, strand, gene_name, cancer, glioblastoma)
setkey(ear045_T3oxBS, chr, start, end, strand, gene_name, cancer, glioblastoma)

M8 <- merge(ear042_M8BS, ear043_M8oxBS, suffixes = c(".BS", ".oxBS"))
T3 <- merge(ear044_T3BS, ear045_T3oxBS, suffixes = c(".BS", ".oxBS"))

setkey(M8, chr, start, end, strand, gene_name, cancer, glioblastoma)
setkey(T3, chr, start, end, strand, gene_name, cancer, glioblastoma)
M8.T3 <- merge(M8, T3, suffixes = c(".M8", ".T3"))

xdata<-M8.T3[gene_name == "ATRX"]

tab<-xdata[, list(pct_5mC.M8=100*sum(cnt_met.oxBS.M8)/sum(cnt_tot.oxBS.M8), pct_5hmC.M8=100*(sum(cnt_met.BS.M8)/sum(cnt_tot.BS.M8) - sum(cnt_met.oxBS.M8)/sum(cnt_tot.oxBS.M8)), pct_C_other.M8=100*(1-sum(cnt_met.BS.M8)/sum(cnt_tot.BS.M8)), pct_5mC.T3=100*sum(cnt_met.oxBS.T3)/sum(cnt_tot.oxBS.T3), pct_5hmC.T3=100*(sum(cnt_met.BS.T3)/sum(cnt_tot.BS.T3) - sum(cnt_met.oxBS.T3)/sum(cnt_tot.oxBS.T3)), pct_C_other.T3=100*(1-sum(cnt_met.BS.T3)/sum(cnt_tot.BS.T3))), by=gene_name]

tab.summary <- reshape2::melt(tab, id.vars = "gene_name", variable.name = "sample", value.name = "pct_met")
tab.summary[, sample2 := factor(c(rep("Margin",3),rep("Tumour",3)))]
tab.summary[, label := factor(rep(c("5mC", "5hmC", "C"),2), levels = c("C", "5hmC", "5mC"))]
tab.summary[, sample := NULL]

gg <- ggplot(tab.summary, aes(x = sample2, y = pct_met, fill = label)) +
geom_bar(stat = "identity") +
xlab("") +
ylab("% Modification") +
theme_classic() +
theme(legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_text(size=10), axis.title.y = element_text(size=10)) +
scale_fill_manual(values=c("lightgray", "#009E73", "#D55E00"))
ggsave("/figures/ATRX.promoter.BS.oxBS.summary.png", width = 7/2.54, height = 4.75/2.54)
```

<img src="figures/ATRX.promoter.BS.oxBS.summary.png" width="300">
