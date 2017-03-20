## Absence of 5hmC in tumor

We ask the question: At sites where 5hmC in margin is high and it is low in tumor, what replaces 5hmC? Does the missing 5hmC results in "C" or "5mC" in tumor?

Here we don't consider the level of 5hmC in tumor as it is generally low. 

```
cd /mnt/nfs/nas/sblab_data1/berald01/projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/

R
library(scales)
library(ggplot2)

tOxBs<- fread('zcat tumor_posterior.txt.gz')
mOxBs<- fread('zcat margin_posterior.txt.gz')
hmcMarg<- fread('zcat hmC_margin_posterior.txt.gz') 

dat<- merge(hmcMarg[, list(chrom, start, mode)], mOxBs[library_id == 'ear043_M8oxBS', list(chrom, start, mode)], by= c('chrom', 'start'), suffixes= c('.hmC_M', '.mC_M'))
dat<- merge(dat, tOxBs[library_id == 'ear045_T3oxBS', list(chrom, start, mode)], by= c('chrom', 'start'))
setnames(dat, 'mode', 'mode.mC_T')

dat[, MmCTier := as.factor(round(mode.mC_M * 10) / 10)]
hline<- data.table(MmCTier= unique(dat$MmCTier), level= as.numeric(as.character(unique(dat$MmCTier))))

gg<- ggplot(data= dat[seq(1, nrow(dat), length.out= 100000)][MmCTier != 1], aes(x= ifelse(mode.hmC_M < 0, 0, ifelse(mode.hmC_M > 0.4, 0.4, mode.hmC_M)), y= mode.mC_T)) +
    geom_point(size= 0.5, alpha= 0.2) +
    geom_smooth() + 
    facet_wrap(~MmCTier, nrow= 3) +
    xlab('% 5hmC Margin') +
    ylab('% 5mC in Tumor') +
    ggtitle('5hmC in margin is 5mC in tumor') +
    geom_hline(data = hline[MmCTier != 1], aes(yintercept = level), colour= 'red', linetype= 'dashed')
ggsave('5hmC_margin_vs_5mC_tumor.pdf', w= 24/2.54, h= 18/2.54)
system('rsync --remove-source-files 5hmC_margin_vs_5mC_tumor.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/')

keep<- seq(1, nrow(dat), length.out= 10000)
gg<- ggplot(data= dat[keep], aes(x= mode.mC_M, y= mode.mC_T, colour= ifelse(mode.hmC_M < 0, 0, ifelse(mode.hmC_M > 0.4, 0.4, mode.hmC_M)))) +
    geom_point(size= 0.7) +
    scale_colour_gradient2('Margin % 5hmC', low=muted("red"), high=muted("blue"), midpoint= 0.2, mid= 'white') +
    geom_abline(intercept= 0, slope= 1, colour= 'red', linetype= 'dashed') +
    xlab('Margin %5mC') +
    ylab('Tumor %5mC') +
    ggtitle('5mC in margin vs tumor')
ggsave('5mC_margin_vs_5mC_tumor.pdf', w= 16/2.54, h= 12/2.54)
system('rsync --remove-source-files 5mC_margin_vs_5mC_tumor.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/')
```