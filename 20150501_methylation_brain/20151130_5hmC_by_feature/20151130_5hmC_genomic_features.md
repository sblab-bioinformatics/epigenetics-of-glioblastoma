## Distribution of 5hmC in margin by genomic feature

We want to see where we find 5hmC in margin (tumour doesn't have much 5hmC, if at all). We use the estimated from 5hmC from bayes method

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/

zcat hmC_margin_posterior.txt.gz \
| awk -v OFS="\t" 'NR > 1 {print $9, $10, $11, $8}' > /tmp/5hmC_Margin.met_diff.bedGraph

zcat hmC_tumor_posterior.txt.gz \
| awk -v OFS="\t" 'NR > 1 {print $9, $10, $11, $8}' > /tmp/5hmC_Tumor.met_diff.bedGraph

for bed in hg19.cpgIslandExt.bed \
           hg19.cpgIslandExtShelves.bed \
           hg19.cpgIslandExtShores.bed \
           hg19.intergenic.bed \
           hg19.refgene.3utr.bed \
           hg19.refgene.5utr.bed \
           hg19.refgene.exons.bed \
           hg19.refgene.introns.bed \
           hg19.refgene.promoters_1000upstream_0downstream.bed
do
echo $bed
bname=`basename $bed .bed`
intersectBed -sorted -a <(sortBed -i ../../reference_data/$bed | mergeBed) -b /tmp/5hmC_Margin.met_diff.bedGraph -wa -wb \
| groupBy -g 1,2,3 -c 1,7 -o count,mean > $bname.margin.tmp.bed

intersectBed -sorted -a <(sortBed -i ../../reference_data/$bed | mergeBed) -b /tmp/5hmC_Tumor.met_diff.bedGraph -wa -wb \
| groupBy -g 1,2,3 -c 1,7 -o count,mean > $bname.tumor.tmp.bed
done

echo "chrom start end n_cpg pct_met feature tissue" | tr ' ' '\t' > /tmp/5hmC_MarginTumor.bed
tableCat.py -i hg19.*.tmp.bed -r '.tmp.bed' | sed 's/hg19\.//' | sed 's/refgene.//' | awk -v OFS="\t" '{gsub(/\./, "\t", $6); print $0}' >> /tmp/5hmC_MarginTumor.bed

rm hg19.*.tmp.bed
rm /tmp/5hmC*.met_diff.bedGraph

R
library(ggplot2)
library(data.table)
library(reshape2)
bed<- fread('/tmp/5hmC_MarginTumor.bed')
bed[, pct_met := 100 * pct_met]
bed[, feature := sub('cpgIslandExtShelves', 'CpGi Shelves', feature)]
bed[, feature := sub('cpgIslandExtShores', 'CpGi Shores', feature)]
bed[, feature := sub('cpgIslandExt', 'CpG islands', feature)]
bed[, feature := sub('promoters_1000upstream_0downstream', 'Promoters', feature)]

dat<- bed
# dat<- bed[feature %in% c("CpG islands", "5utr", "3utr", "exons", "introns", "intergenic", "Promoters")]

## Ordering of boxes by oxBS median
## --------------------------------
mdn<- dat[tissue == 'margin', list(median= median(pct_met)), by= list(feature)][order(median)]
dat[, feature := factor(dat$feature, levels= mdn$feature)]

gg<- ggplot(data= dat, aes(x= feature, y= pct_met, by= tissue, colour= tissue)) +
   geom_boxplot(scale= 'area', position=position_dodge(0.75), outlier.size= 0.75, outlier.colour= 'grey30') +
   geom_hline(yintercept= 0, colour= 'blue', alpha= 0.2, size= 1) + 
   ggtitle("5hmC by genomic feature") + 
   ylab('% 5hmC') +
   ylim(-30, 60) + 
   theme(axis.text.x= element_text(angle= 45, hjust= 1))
ggsave('boxplot_pct_5hmC_chr18_cpg_by_region_TumorMargin.pdf', w= 18, h= 12, units= 'cm')
system(
"rsync --remove-source-files boxplot_pct_5hmC_chr18_cpg_by_region_TumorMargin.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/")

quit(save= 'no')
rm /tmp/5hmC_MarginTumor.bed
```
