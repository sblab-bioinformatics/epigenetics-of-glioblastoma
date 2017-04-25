
This script is an analysis to respond to Shankar:

In this dataset, whilst we saw global hmC loss (also reflected across all region of gene structure) in the tumour.  Were there any genes for which hmC levels actually went up in the tumour relative to margin? I ask because, the literature suggests that in transitional cases (e.g. development) that tissue-specific (and gene-upregulation specific) hmC levels change and can become higher in gene bodies.  For cancers, the dogma is global hmC loss, but if the transcriptional hypothesis (i.e. hmC in gene bodies marks transcriptionally active genes) is correct, there should be some cases of hmC elevation in gene bodies in the tumour. As a guide, it might be worth considering cancer genes unregulated in the tumour, to narrow the list of genes to look at.

Our analysis in the paper was restricted to promoters.  The literature appears to be stronger for gene body hmC levels vs transcription.

It is all based on our previous work:

- https://github.com/sblab-bioinformatics/epigenetics-of-glioblastoma/blob/master/20150501_methylation_brain/20160226_tsg_ocg.md
- https://github.com/sblab-bioinformatics/epigenetics-of-glioblastoma/blob/master/20150501_methylation_brain/20160531_figure_editing.md#analysis-of-expression-of-top-10-5hmc-differences-between-tumour-and-margin



## Gene body definition

Gene body was defined as the entire gene from the transcription start site to the end of the transcript [here](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018844)

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap
mkdir genebody_5hmC
cd genebody_5hmC
```

https://www.biostars.org/p/103411/

```r
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/lustre/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf", format="gtf")
genes <- genes(txdb)
write.table(as.data.frame(genes)[,-4], file = "/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/genebody_5hmC/20170425_hg19_genecoords.bed", row.names = F, col.names = F, sep = "\t", quote = F)

```

After searching TET in IGV, the gene boundaries look ok.

Sort:

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/genebody_5hmC
sort -k1,1 -k2,2n 20170425_hg19_genecoords.bed > 20170425_hg19_genecoords_sorted.bed
wc -l 20170425_hg19_genecoords_sorted.bed # 22833
rm 20170425_hg19_genecoords.bed
```



## Intersect gene bodies and 5mC/5hmC levels

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/genebody_5hmC

bedtools intersect -a ../methylation_cpg/20160210/ear042_M8BS.cpg.bedGraph -b 20170425_hg19_genecoords_sorted.bed -wa -wb -sorted | cut -f 1-6,10-11 > 20170425_ear042_M8BS.cpg_hg19_genecoords_sorted.bed
bedtools intersect -a ../methylation_cpg/20160210/ear043_M8oxBS.cpg.bedGraph -b 20170425_hg19_genecoords_sorted.bed -wa -wb -sorted | cut -f 1-6,10-11 > 20170425_ear043_M8oxBS.cpg_hg19_genecoords_sorted.bed
bedtools intersect -a ../methylation_cpg/20160210/ear044_T3BS.cpg.bedGraph -b 20170425_hg19_genecoords_sorted.bed -wa -wb -sorted | cut -f 1-6,10-11 > 20170425_ear044_T3BS.cpg_hg19_genecoords_sorted.bed
bedtools intersect -a ../methylation_cpg/20160210/ear045_T3oxBS.cpg.bedGraph -b 20170425_hg19_genecoords_sorted.bed -wa -wb -sorted | cut -f 1-6,10-11 > 20170425_ear045_T3oxBS.cpg_hg19_genecoords_sorted.bed

wc -l 20170425_ear*.bed
#  28082863 20170425_ear042_M8BS.cpg_hg19_genecoords_sorted.bed
#  28039066 20170425_ear043_M8oxBS.cpg_hg19_genecoords_sorted.bed
#  28061617 20170425_ear044_T3BS.cpg_hg19_genecoords_sorted.bed
#  28052114 20170425_ear045_T3oxBS.cpg_hg19_genecoords_sorted.bed

wc -l ../methylation_cpg/20160210/ear*.cpg.bedGraph
#  54834640 ../methylation_cpg/20160210/ear042_M8BS.cpg.bedGraph
#  54765875 ../methylation_cpg/20160210/ear043_M8oxBS.cpg.bedGraph
#  54801746 ../methylation_cpg/20160210/ear044_T3BS.cpg.bedGraph
#  54788215 ../methylation_cpg/20160210/ear045_T3oxBS.cpg.bedGraph

cut -f 8 20170425_ear042_M8BS.cpg_hg19_genecoords_sorted.bed | sort | uniq | wc -l # 22338 genes in total
cut -f 8 20170425_ear043_M8oxBS.cpg_hg19_genecoords_sorted.bed | sort | uniq | wc -l # 22345 genes in total
cut -f 8 20170425_ear044_T3BS.cpg_hg19_genecoords_sorted.bed | sort | uniq | wc -l # 22346 genes in total
cut -f 8 20170425_ear045_T3oxBS.cpg_hg19_genecoords_sorted.bed | sort | uniq | wc -l # 22345 genes in total

```

The columns of the files generated above are:

- 1 chr
- 2 start
- 3 end
- 4 cnt_met
- 5 cnt_tot
- 6 + or - (C or G)
- 7 + or - (gene strandness)
- 8 gene name



## Transcript levels vs. 5mC/5hmC levels in gene bodies

https://github.com/sblab-bioinformatics/epigenetics-of-glioblastoma/blob/master/20150501_methylation_brain/20160531_figure_editing.md#analysis-of-expression-of-top-10-5hmc-differences-between-tumour-and-margin
