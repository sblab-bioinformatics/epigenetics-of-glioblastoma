

The goal of this script is to normalise the figures and panels that are going to go inside the manuscript.

### 5mC distribution by genomic features

Modified from:

https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160126_methyl_distr/20160126_5mC_distr.md

Copy input files to lustre. In nas-srv001:

```bash
cp /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160126_methyl_distr/* /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/data

```


```R
library(ggplot2)
library(data.table)

bed <- fread('zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/data/5mC_MarginTumor.bed.gz')
bed$feature <- factor(bed$feature)
levels(bed$feature)
levels(bed$feature) <- c("3'-UTR", "5'-UTR", "CpGi", "CpGi shelves", "CpGi shores", "Exons", "Intergenic", "Introns", "Promoters")
table(bed$feature)
bed$tissue <- factor(bed$tissue)
levels(bed$tissue)
levels(bed$tissue) <- c("Margin", "Tumour")
table(bed$tissue)

dat <- bed
mdn<- dat[tissue == 'Margin', list(median = median(pct_met)), by = list(feature)][order(median)]
dat[, feature := factor(dat$feature, levels= mdn$feature)]

gg <- ggplot(data = dat, aes(x = feature, y = pct_met, fill = tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75) +
theme_classic() +
xlab("") +
ylab('% 5mC') +
scale_fill_manual(values=c("green4", "red3")) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank())
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160531_5mC_region.pdf', width = 14/2.54, height = 10/2.54)

```



### 5hmC distribution by genomic features


```R
library(ggplot2)
library(data.table)

bed <- fread('zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/data/5hmC_MarginTumor.bed.gz')
bed$feature <- factor(bed$feature)
levels(bed$feature)
levels(bed$feature) <- c("3'-UTR", "5'-UTR", "CpGi", "CpGi shelves", "CpGi shores", "Exons", "Intergenic", "Introns", "Promoters")
table(bed$feature)
bed$tissue <- factor(bed$tissue)
levels(bed$tissue)
levels(bed$tissue) <- c("Margin", "Tumour")
table(bed$tissue)

dat <- bed
dat[, feature := factor(dat$feature, levels= c("CpGi", "Promoters", "5'-UTR", "CpGi shores", "3'-UTR", "CpGi shelves", "Exons", "Introns", "Intergenic"))]

gg <- ggplot(data = dat, aes(x = feature, y = pct_met, fill = tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75) +
theme_classic() +
xlab("") +
ylab('% 5hmC') +
scale_fill_manual(values=c("green4", "red3")) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank())
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160622_5hmC_region.pdf', width = 14/2.54, height = 10/2.54)

```



### Merge 5mC and 5hmC distributions by genomic features in one single plot

```R
library(ggplot2)
library(data.table)

bed_5mC <- fread('zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/data/5mC_MarginTumor.bed.gz')
bed_5mC$feature <- factor(bed_5mC$feature)
levels(bed_5mC$feature)
levels(bed_5mC$feature) <- c("3'-UTR", "5'-UTR", "CpGi", "CpGi shelves", "CpGi shores", "Exons", "Intergenic", "Introns", "Promoters")
table(bed_5mC$feature)
bed_5mC$tissue <- factor(bed_5mC$tissue)
levels(bed_5mC$tissue)
levels(bed_5mC$tissue) <- c("Margin", "Tumour")
table(bed_5mC$tissue)

dat_5mC <- bed_5mC
mdn_5mC<- dat_5mC[tissue == 'Margin', list(median = median(pct_met)), by = list(feature)][order(median)]
dat_5mC[, feature := factor(dat_5mC$feature, levels= mdn_5mC$feature)]

bed_5hmC <- fread('zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/data/5hmC_MarginTumor.bed.gz')
bed_5hmC$feature <- factor(bed_5hmC$feature)
levels(bed_5hmC$feature)
levels(bed_5hmC$feature) <- c("3'-UTR", "5'-UTR", "CpGi", "CpGi shelves", "CpGi shores", "Exons", "Intergenic", "Introns", "Promoters")
table(bed_5hmC$feature)
bed_5hmC$tissue <- factor(bed_5hmC$tissue)
levels(bed_5hmC$tissue)
levels(bed_5hmC$tissue) <- c("Margin", "Tumour")
table(bed_5hmC$tissue)

dat_5hmC <- bed_5hmC
dat_5hmC[, feature := factor(dat_5hmC$feature, levels= c("CpGi", "Promoters", "5'-UTR", "CpGi shores", "3'-UTR", "CpGi shelves", "Exons", "Introns", "Intergenic"))]

l = list(dat_5mC, dat_5hmC)
setattr(l, 'names', c("5mC", "5hmC"))
dat_5mC_5hmC <- rbindlist(l, use.names=TRUE, idcol="ID")

gg <- ggplot(data = dat_5mC_5hmC, aes(x = feature, y = pct_met, fill = tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75, lwd = 0.25, aes(colour = factor(ID, levels(factor(ID))[c(2,1)]))) +
theme_classic() +
xlab("") +
ylab('% Modification') +
scale_fill_manual(values=c("green4", "red3")) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
scale_colour_manual(values = c("turquoise3", "green"))
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160701_5mC_5hmC_region.pdf', width = 16/2.54, height = 12/2.54)


# 20160722 update
dat_5mC_5hmC[, modification_tissue := paste(ID, tissue, sep = " ")]
table(dat_5mC_5hmC$modification_tissue)

gg <- ggplot(data = dat_5mC_5hmC, aes(x = feature, y = pct_met, fill = modification_tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75, lwd = 0.25) +
theme_classic() +
xlab("") +
ylab('% Modification') +
scale_fill_manual(values=c("green", "green4", "turquoise3", "turquoise4")) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank())
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160722_5mC_5hmC_region.pdf', width = 16/2.54, height = 12/2.54)


# 20160823 update
dat_5mC_5hmC[, modification_tissue := paste(ID, tissue, sep = " ")]
table(dat_5mC_5hmC$modification_tissue)
dat_5mC_5hmC[, modification_tissue := factor(modification_tissue, levels = c("5mC Tumour", "5mC Margin", "5hmC Margin", "5hmC Tumour"))]

# [Wong2011](http://www.nature.com/nmeth/journal/v8/n6/full/nmeth.1618.html)
# 5mC Margin Vermillion rgb(213, 94, 0, max = 255) "#D55E00"
# 5mC Tumour Vermillion (manually darkened using Inkscape) "#A24800" rgb(162, 72, 0, max = 255)
# 5hmC Margin Bluish green rgb(0, 158, 115, max = 255) "#009E73"
# 5hmC Tumour Bluish green (manually darkened using Inkscape) "#006C4F" rgb(0, 108, 79, max = 255)

gg <- ggplot(data = dat_5mC_5hmC, aes(x = feature, y = pct_met, fill = modification_tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75, lwd = 0.25) +
theme_classic() +
xlab("") +
ylab('% Modification') +
scale_fill_manual(values=c("#A24800", "#D55E00", "#009E73", "#006C4F")) +
theme(legend.title = element_blank()) +
scale_x_discrete(labels = c("CpGi", "Promoters", "5'-UTR", "CpGi\nshores", "3'-UTR", "CpGi\nshelves", "Exons", "Introns", "Intergenic"))
#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) # vertical labels x ticks
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160823_5mC_5hmC_region.pdf', width = 20/2.54, height = 12/2.54)


# 20160902 update
dat_5mC_5hmC[, modification_tissue := paste(ID, tissue, sep = " ")]
table(dat_5mC_5hmC$modification_tissue)
dat_5mC_5hmC[, modification_tissue := factor(modification_tissue, levels = c("5mC Tumour", "5mC Margin", "5hmC Margin", "5hmC Tumour"))]

# [Wong2011](http://www.nature.com/nmeth/journal/v8/n6/full/nmeth.1618.html)
# 5mC Margin Vermillion rgb(213, 94, 0, max = 255) "#D55E00"
# 5mC Tumour Vermillion (manually darkened using Inkscape) "#A24800" rgb(162, 72, 0, max = 255)
# 5hmC Margin Bluish green rgb(0, 158, 115, max = 255) "#009E73"
# 5hmC Tumour Bluish green (manually darkened using Inkscape) "#006C4F" rgb(0, 108, 79, max = 255)

gg <- ggplot(data = dat_5mC_5hmC, aes(x = feature, y = pct_met, fill = modification_tissue)) +
geom_boxplot(outlier.shape = NA, position=position_dodge(0.75), width = 0.75, lwd = 0.25) +
theme_classic() +
xlab("") +
ylab('% Modification') +
scale_fill_manual(values=c("#A24800", "#D55E00", "#009E73", "#006C4F")) +
theme(legend.title = element_blank(), legend.position = "top", axis.text.x = element_text(size=8), axis.title.y = element_text(size=10)) +
scale_x_discrete(labels = c("CpGi", "Promoters", "5'-UTR", "CpGi\nshores", "3'-UTR", "CpGi\nshelves", "Exons", "Introns", "Intergenic"))
ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160902_5mC_5hmC_region.pdf', width = 13.5/2.54, height = 12/2.54)

```




### Absence of 5hmC in tumour

Modified from:

https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151102_Tumor_vs_Margin_by_bayesian_approach/20151210_missing_hmC_in_tumor.md

We could do it just reproducing the plots using the bayesian approach by Dario's. However this is chr18 only and it will be better to have a sample of CpG sites from all the genome. Therefore we could sample sites from the substraction and plot these only

Copy input files to lustre. In nas-srv001:

```R
library(data.table)
library(ggplot2)
library(scales)

ear042_M8BS <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear042_M8BS.cpg.bedGraph")
setnames(ear042_M8BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear043_M8oxBS <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear043_M8oxBS.cpg.bedGraph")
setnames(ear043_M8oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear044_T3BS <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear044_T3BS.cpg.bedGraph")
setnames(ear044_T3BS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

ear045_T3oxBS <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear045_T3oxBS.cpg.bedGraph")
setnames(ear045_T3oxBS, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))

setkey(ear042_M8BS, chr, start, end, strand)
setkey(ear043_M8oxBS, chr, start, end, strand)
setkey(ear044_T3BS, chr, start, end, strand)
setkey(ear045_T3oxBS, chr, start, end, strand)

M8 <- merge(ear042_M8BS, ear043_M8oxBS, suffixes = c(".BS", ".oxBS"))
T3 <- merge(ear044_T3BS, ear045_T3oxBS, suffixes = c(".BS", ".oxBS"))
rm(ear042_M8BS, ear043_M8oxBS, ear044_T3BS, ear045_T3oxBS)

setkey(M8, chr, start, end, strand)
setkey(T3, chr, start, end, strand)
M8.T3 <- merge(M8, T3, suffixes = c(".M8", ".T3"))
rm(M8, T3)

M8.T3[, pct_5mC.M8 := 100*cnt_met.oxBS.M8/cnt_tot.oxBS.M8][, pct_5hmC.M8 := 100*(cnt_met.BS.M8/cnt_tot.BS.M8 - cnt_met.oxBS.M8/cnt_tot.oxBS.M8)][, pct_5mC.T3 := 100*cnt_met.oxBS.T3/cnt_tot.oxBS.T3][, pct_5hmC.T3 := 100*(cnt_met.BS.T3/cnt_tot.BS.T3 - cnt_met.oxBS.T3/cnt_tot.oxBS.T3)]

par(mfrow=c(2,1))
plot(density(M8.T3$cnt_tot.oxBS.M8), xlim = c(0,100))
plot(density(M8.T3$cnt_tot.oxBS.T3), xlim = c(0,100))

# clearly need to sample from sites with more than 20 counts
M8.T3.20 <- M8.T3[cnt_tot.oxBS.M8 > 20 & cnt_tot.oxBS.T3 > 20]

set.seed(1)
M8.T3.20.sample <- M8.T3.20[sample(c(1:nrow(M8.T3.20)), 100000)]


# Scatterplot
gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC\nMargin", midpoint = 25, low = muted("red"), high = muted("blue"), mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'blue', linetype= 'dashed')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160609_5mC.M8.5mC.T3.5hmC.M8.scatterplot.pdf", width = 16/2.54, height = 12/2.54)
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160609_5mC.M8.5mC.T3.5hmC.M8.scatterplot.png", width = 16/2.54, height = 12/2.54)


# Boxplot
mean(M8.T3.20.sample$pct_5mC.M8) # 53.49721
mean(M8.T3.20.sample$pct_5mC.T3) # 65.06919
M8.T3.20.sample[, label.5mC := ifelse(pct_5mC.M8 > mean(pct_5mC.M8), ifelse(pct_5mC.T3 > mean(pct_5mC.T3), "h5mC_h5mC", "h5mC_l5mC"), ifelse(pct_5mC.T3 > mean(pct_5mC.T3), "l5mC_h5mC", "l5mC_l5mC"))]
table(M8.T3.20.sample$label.5mC)
#h5mC_h5mC h5mC_l5mC l5mC_h5mC l5mC_l5mC
#    51288      8085     15459     25168

gg2 <- ggplot(M8.T3.20.sample, aes(x = factor(label.5mC, levels(factor(label.5mC))[c(3,4,1,2)]), y = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) + geom_boxplot(outlier.shape = NA) + theme_classic() + xlab("") + ylab("% 5hmC Margin") + scale_x_discrete(labels = c("low\nhigh", "low\nlow", "high\nhigh", "high\nlow"))
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160609_5mC.M8.5mC.T3.5hmC.M8.boxplot.pdf", width = 10/2.54, height = 10/2.54)


# 20160823 update

# Scatterplot
gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC Margin", midpoint = 25, low = "#D55E00", high = "#009E73", mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dashed')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160823_5mC.M8.5mC.T3.5hmC.M8.scatterplot.pdf", width = 16/2.54, height = 12/2.54)

gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC Margin", midpoint = 25, low = "#D55E00", high = "#009E73", mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dashed') +
theme(legend.position="none")
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160823_5mC.M8.5mC.T3.5hmC.M8.scatterplot.png", width = 12/2.54, height = 12/2.54)



# 20160824 update - change to a greenscale colour palette

# Scatterplot
gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC Margin", midpoint = 5, low = "white", high = "#009E73", mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dashed')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160824_5mC.M8.5mC.T3.5hmC.M8.scatterplot.greenscale.pdf", width = 16/2.54, height = 12/2.54)

gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC Margin", midpoint = 5, low = "white", high = "#009E73", mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dashed') +
theme(legend.position="none")
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160824_5mC.M8.5mC.T3.5hmC.M8.scatterplot.greenscale.png", width = 12/2.54, height = 12/2.54)


# Scatterplot - just playing with scale_colour_gradientn
gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradientn("% 5hmC Margin", colours = colorRampPalette(c("white", "white", "#009E73"), bias = 3)(100)) +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dashed')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160824_5mC.M8.5mC.T3.5hmC.M8.scatterplot.greenscale.bias3.pdf", width = 16/2.54, height = 12/2.54)



# 20160902 update - keep green and orange scale instead

# Scatterplot
gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
theme(axis.title = element_text(size=10), axis.text = element_text(size=10), legend.position = "top",  legend.title = element_text(size=10)) +
scale_colour_gradient2("% 5hmC Margin", midpoint = 25, low = "#D55E00", high = "#009E73", mid = 'white') +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dotted')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160902_5mC.M8.5mC.T3.5hmC.M8.scatterplot.pdf", width = 11/2.54, height = 12/2.54)

gg <- ggplot(M8.T3.20.sample, aes(x = pct_5mC.M8, y = pct_5mC.T3, colour = ifelse(pct_5hmC.M8 < 0, 0, ifelse(pct_5hmC.M8 > 50, 50, pct_5hmC.M8)))) +
geom_point(size=0.05) +
xlab("% 5mC Margin") +
ylab("% 5mC Tumour") +
theme_classic() +
scale_colour_gradient2("% 5hmC Margin", midpoint = 25, low = "#D55E00", high = "#009E73", mid = 'white') +
theme(axis.title = element_text(size=10), axis.text = element_text(size=10),legend.position="none") +
geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), colour= 'black', linetype= 'dotted')
ggsave("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160902_5mC.M8.5mC.T3.5hmC.M8.scatterplot.png", width = 11/2.54, height = 10/2.54)


```




### Analysis of expression of top 10 5hmC differences between tumour and margin

I followed Dario's code (also introduced a few small changes) to find out which are the gene promoters with top 10 5hmC differences between tumour and margin
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160705_top_methyl_proms/20160705_top_methyl_proms.md


```R
library(data.table)
library(ggplot2)


## Promoter 5(h)mC
proms<- fread('tableCat.py -r "\\..*" -i /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/tsgenes_oncogenes/data/methylation_cpg_promoter/ear*.cpg_genes.promoters.1000.sorted.bed')
setnames(proms, names(proms), c('chrom', 'start', 'end', 'cnt_met', 'cnt_tot', 'strand', 'gene_name', 'unused', 'library_id'))
proms[, unused := NULL]
proms<- proms[cnt_tot > 0,]
proms<- proms[, list(cnt_met= sum(cnt_met), cnt_tot= sum(cnt_tot), .N), by= list(gene_name, library_id)]


## ------[ Promoters to include in further analysis ] -----
## Number of cytosines per promoter. MEMO: These are not CpG
ggplot(data= proms, aes(x= ifelse(N > 100, 100, N))) + geom_histogram() + facet_wrap(~library_id)
## Promoter depth
ggplot(data= proms, aes(x= ifelse(cnt_tot > 1000, 1000, cnt_tot))) + geom_histogram() + facet_wrap(~library_id)

# ncyt gives the number of libraries with more than x cytosines in the promoter.
# ndepth is the num of libs with cnt_tot > x
promset<- proms[, list(ncyt= sum(N > 20), ndepth= sum(cnt_tot > 100)), by= list(gene_name)]
stopifnot(nrow(promset) == length(unique(proms$gene_name)))


## Consider only these promoters:
nlibs<- length(unique(proms$library_id))
promLst<- promset[ncyt == nlibs & ndepth == nlibs, gene_name]


## Prepare table of 5mC and 5hmC in tumor and margin
proms[, pct_met := 100 * (cnt_met/cnt_tot)]
prommet<- dcast.data.table(data= proms[gene_name %in% promLst], gene_name ~ library_id, value.var= 'pct_met')
prommet<- prommet[, list(gene_name,
               mar_mc= ear043_M8oxBS,  
               mar_hmc= ear042_M8BS - ear043_M8oxBS,
               tum_mc= ear044_T3BS,
               tum_hmc= ear044_T3BS - ear045_T3oxBS)]
prommet[, mar_hmc := ifelse(mar_hmc < 0, 0, mar_hmc)]
prommet[, tum_hmc := ifelse(tum_hmc < 0, 0, tum_hmc)]


## Table Gene Expression
tx<- fread('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/rnaseq/data/tx_quant_lfc.bed')

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


## Few gene names appear on forward and reverse strand! Get rid of them:
dups<- gx[duplicated(gx[, Associated_Gene_Name]), Associated_Gene_Name] ## "KBTBD4" "NPIPA7" "ZNF668" "LIMS3" "RBL1" "CKS1B"
gx<- gx[!Associated_Gene_Name %in% dups]


## 5(h)mC vs Expr
metexpr<- merge(prommet, gx, by.x= 'gene_name', by.y= 'Associated_Gene_Name')
stopifnot(nrow(metexpr) == length(unique(metexpr$gene_name)))

gg<- ggplot(data= metexpr, aes(x= tum_mc - mar_mc, y= tpmLog2FC,
        col= densCols(tum_mc - mar_mc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
        )
    ) +
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    theme_classic() +
    xlab("% 5mC promoter Tumour - % 5mC promoter Margin") +
    ylab(expression("log"[2]*"FC"))

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160719_promoter_5mC_difference_log2FC.pdf', width= 14, height= 12, units= 'cm')

gg<- ggplot(data= metexpr, aes(x= tum_hmc - mar_hmc, y= tpmLog2FC,
        col= densCols(metexpr$tum_hmc - metexpr$mar_hmc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))) +
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    theme_classic() +
    xlab("% 5hmC promoter Tumour - % 5hmC promoter Margin") +
    ylab(expression("log"[2]*"FC"))

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160719_promoter_5hmC_difference_log2FC.pdf', width= 14, height= 12, units= 'cm')

gg<- ggplot(data= metexpr, aes(x= (tum_hmc/tum_mc) - (mar_hmc/mar_mc), y= tpmLog2FC,
    col= densCols((tum_hmc/tum_mc) - (mar_hmc/mar_mc), tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))) +
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    theme_classic() +
    xlab("(% 5hmC/% 5mC) promoter Tumour - (% 5hmC/% 5mC) promoter Margin") +
    ylab(expression("log"[2]*"FC"))

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160719_promoter_5hmC_5mC_ratio_difference_log2FC.pdf', width= 14, height= 12, units= 'cm')


## Annotate promoters
library(biomaRt)
# Hosts can be found at http://www.ensembl.org/info/website/archives/index.html
mart<- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "dec2015.archive.ensembl.org", dataset= 'hsapiens_gene_ensembl')
genedescr<- data.table(getBM(
    filters= "external_gene_name",
    attributes= c("external_gene_name", "description"),
    values= metexpr$gene_name,
    mart=mart))
genedescr<- genedescr[gsub("^//s+|//s+$", "", description) != ""] ## Remove genes w/o descr
genedescr<- genedescr[, list(description= paste(description, collapse= "|")), by= external_gene_name] ## Collapse dup descr
stopifnot(length(genedescr$external_gene_name) == length(unique(genedescr$external_gene_name)))

xmrg<- merge(metexpr, genedescr, by.x= 'gene_name', by.y= 'external_gene_name', all.x= TRUE)


## Top 10 gene promoters with larger 5hmC differences between tumour and margin
xmrg[, diff_hmc := mar_hmc - tum_hmc]
xmrg[order(-diff_hmc)][c(1:10)]
# PIK3C2B might be of interest
xmrg[order(diff_hmc)][c(1:10)]


## 20160722 - add names to top genes
gg<- ggplot(data= metexpr, aes(x= mar_hmc - tum_hmc, y= tpmLog2FC,
        col= densCols(metexpr$mar_hmc - metexpr$tum_hmc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))), label = gene_name)) +
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    theme_classic() +
    xlab("% 5hmC promoter Margin - % 5hmC promoter Tumour") +
    ylab(expression("log"[2]*"FC")) +
    geom_text(aes(label = ifelse((mar_hmc - tum_hmc) > 37, as.character(gene_name), ifelse(((mar_hmc - tum_hmc) > 30) & ((tpmLog2FC > 5) | (tpmLog2FC < -5)), as.character(gene_name), ""))), hjust=-0.1, vjust=-0.1, size=3)

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160722_promoter_5hmC_difference_log2FC_genenames.pdf', width= 14, height= 12, units= 'cm')


## with ggrepel (https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
# install.packages("/home/martin03/tmp/building/ggrepel_0.5.tar.gz", repos = NULL, type="source")
library(ggrepel)
metexpr[, diff_hmc := mar_hmc - tum_hmc]
gg<- ggplot(data= metexpr, aes(x= mar_hmc - tum_hmc, y= tpmLog2FC)) +
    geom_point(size= 0.1, color = densCols(metexpr$mar_hmc - metexpr$tum_hmc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    theme_classic() +
    xlab("% 5hmC promoter Margin - % 5hmC promoter Tumour") +
    ylab(expression("log"[2]*"FC")) +
    geom_label_repel(data = metexpr[order(-diff_hmc)][c(1:10)], aes(label = gene_name), size=3, force = 1)

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160722_promoter_5hmC_difference_log2FC_genenames_ggrepel.pdf', width= 14, height= 12, units= 'cm')


## 20160907 - repeat ggrepel plot but adjusting it to Fig. 2e
gg<- ggplot(data= metexpr, aes(x= diff_hmc, y= tpmLog2FC)) +
geom_point(size= 0.1, color = densCols(metexpr$diff_hmc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) +
geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
theme_classic() +
xlab("% 5hmC promoter Margin - % 5hmC promoter Tumour") +
ylab(expression("log"[2]*"FC")) +
geom_label_repel(data = metexpr[order(-diff_hmc)][c(1:10)], aes(label = gene_name), size=3, force = 1) +
theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11), axis.title.x = element_text(size=9.5), axis.title.y = element_text(size=13))

ggsave('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/figures/20160907_promoter_5hmC_difference_log2FC_genenames_ggrepel.pdf', width= 10, height= 8.5, units= 'cm')







## Top 100 gene promoters with larger 5hmC differences between tumour and margin
## Write to table for Euni to digest
xmrg[, diff_hmc := mar_hmc - tum_hmc]
write.table(xmrg[order(-diff_hmc)][c(1:100)], file = "/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/figure_editing/tables/20160722_top100genes_5hmCvsexpression_tumourvsmargin.txt", quote = FALSE, sep = "\t", row.names=FALSE)



## Functional analysis of top genes
## https://www.bioconductor.org/packages/release/BiocViews.html#___GeneSetEnrichment
## limma, CompGo, edgeR, gage, mdgsa, goseq

library(data.table)

proms<- fread('tableCat.py -r "\\..*" -i /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/tsgenes_oncogenes/data/methylation_cpg_promoter/ear*.cpg_genes.promoters.1000.sorted.bed')
setnames(proms, names(proms), c('chrom', 'start', 'end', 'cnt_met', 'cnt_tot', 'strand', 'gene_name', 'unused', 'library_id'))
proms[, unused := NULL]
proms<- proms[cnt_tot > 0,]
proms<- proms[, list(cnt_met= sum(cnt_met), cnt_tot= sum(cnt_tot), .N), by= list(gene_name, library_id)]

promset<- proms[, list(ncyt= sum(N > 20), ndepth= sum(cnt_tot > 100)), by= list(gene_name)]
stopifnot(nrow(promset) == length(unique(proms$gene_name)))

nlibs<- length(unique(proms$library_id))
promLst<- promset[ncyt == nlibs & ndepth == nlibs, gene_name]

proms[, pct_met := 100 * (cnt_met/cnt_tot)]
prommet<- dcast.data.table(data= proms[gene_name %in% promLst], gene_name ~ library_id, value.var= 'pct_met')
prommet<- prommet[, list(gene_name,
               mar_mc= ear043_M8oxBS,  
               mar_hmc= ear042_M8BS - ear043_M8oxBS,
               tum_mc= ear044_T3BS,
               tum_hmc= ear044_T3BS - ear045_T3oxBS)]
prommet[, mar_hmc := ifelse(mar_hmc < 0, 0, mar_hmc)]
prommet[, tum_hmc := ifelse(tum_hmc < 0, 0, tum_hmc)]

tx<- fread('/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/rnaseq/data/tx_quant_lfc.bed')

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

dups<- gx[duplicated(gx[, Associated_Gene_Name]), Associated_Gene_Name] ## "KBTBD4" "NPIPA7" "ZNF668" "LIMS3" "RBL1" "CKS1B"
gx<- gx[!Associated_Gene_Name %in% dups]

metexpr<- merge(prommet, gx, by.x= 'gene_name', by.y= 'Associated_Gene_Name')
stopifnot(nrow(metexpr) == length(unique(metexpr$gene_name)))

metexpr[, diff_hmc := mar_hmc - tum_hmc]


## Get gene ids

library(biomaRt)

listMarts(host="www.ensembl.org")

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset= 'hsapiens_gene_ensembl', host = "jul2016.archive.ensembl.org") # http://www.ensembl.org/info/website/archives/index.html
#listAttributes(mart)
#listFilters(mart)
geneid <- data.table(getBM(attributes = c("external_gene_name", "entrezgene"), values = metexpr$gene_name, mart = mart))
geneid <- geneid[!is.na(geneid$entrezgene),]
geneid <- geneid[, list(entrezgene = paste(entrezgene, collapse= "|")), by= external_gene_name]
stopifnot(length(geneid$external_gene_name) == length(unique(geneid$external_gene_name)))

xmrg<- merge(metexpr, geneid, by.x= 'gene_name', by.y= 'external_gene_name', all.x= TRUE)



## edgeR

library(edgeR)

de <- xmrg[order(-diff_hmc)][c(1:100)]


# GO analysis
go <- goana(de = de$entrezgene, universe = xmrg$entrezgene)
topGO(go, ont = "MF", number = 20)
# MF: the top 100 genes showing clear loss of promoter 5hmC in tumour are linked to the transport of anions
topGO(go, ont = "BP", number = 20)
# BP: immune response, negative regulation of Ras protein signal transduction
topGO(go, ont = "CC", number = 20)
# CC: no clear pattern

# Obtain genes related to each relevant GO term, e.g. GO:0015106
library(org.Hs.eg.db)

org.Hs.eg()
#DB schema: HUMAN_DB
#DB schema version: 2.1
#Organism: Homo sapiens
#Date for NCBI data: 2015-Sep27
#Date for GO data: 20150919
#Date for KEGG data: 2011-Mar15
#Date for Golden Path data: 2010-Mar22
#Date for Ensembl data: 2015-Jul16

org.Hs.eg_dbInfo()

x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0015106"
EG <- mappedLkeys(x)
EG

de[de$entrezgene %in% EG,]


# KEGG analysis
kegg <- kegga(de = de$entrezgene, universe = xmrg$entrezgene)
topKEGG(kegg, number = 20)


```







The volcano plot and the glioblastoma genes plot are here:
https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160303_rnaseq/20160303_rnaseq.md#expression-vs-methylation
