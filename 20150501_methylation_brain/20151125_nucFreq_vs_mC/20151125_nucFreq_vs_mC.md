## Nucleotide frequency vs methylation

We want to see if there is a correlation between nucleotide frequencies and methylation at CpG sites. For margin sample we could look at
5hmC as well. Restrict to chr18 for now.

For this we need:

* Table of methylation. See [20151102 Tumor vs Margin by bayesian approach](20151102_Tumor_vs_Margin_by_bayesian_approach)
* Table of nucleotide freqs from Sergio. See [20151124 SNV analysis](20151124_SNV_analysis)

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/

## 5mC in tumor and margin
../20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/margin_posterior.txt.gz
../20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/tumor_posterior.txt.gz

## Nuc freq in tumor and margin
scp uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_C01.18.txt ./ # Margin
scp uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_E01.18.txt ./ # Tumor

## Convert chrom name and make it bed. Positions appear to be 0-based
## Columns are A, C, G, T
sed 's/^18\t/chr18\t/' LP2000729-DNA_C01.18.txt | awk -v OFS="\t" 'NR>1 {print $1, $2, $2+1, $3, $4, $5, $6, $7}' > margin_DNA.chr18.txt
sed 's/^18\t/chr18\t/' LP2000729-DNA_E01.18.txt | awk -v OFS="\t" 'NR>1 {print $1, $2, $2+1, $3, $4, $5, $6, $7}' > tumor_DNA.chr18.txt
rm LP2000729-DNA_*.18.txt

## MEMO: Position of nucs:
# $5= A
# $6= C
# $7= G
# $8= T

for tissue in margin tumor
do
    ## At each C/G count the C/G (i.e. non SNV), the T/A (i.e. the SNV that bias methylation call) and
    ## the non C/G, T/A (SNV)
    awk -v OFS="\t" '{
    if($4 == "C")
        {c=$6;
         t=$8
         n=$5+$7}
    else if($4 == "G")
        {c=$7;
         t=$5;
         n=$6+$8}
    else
        {print 'WRONG'; exit 1};
    if((c+t+n) != $5+$6+$7+$8){
        {print 'WRONG-2'; exit 1};
    }
    if((c+t+n) > 0)
        {print $1, $2, $3, c, t, n}}' ${tissue}_DNA.chr18.txt > ${tissue}_DNA.chr18.pct_alt.bedGraph &&
    
    ## Methylation
    zcat ../20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/${tissue}_posterior.txt.gz \
    | awk -v OFS="\t" 'NR > 1 && $9 ~ /oxBS/ {print $10, $11, $12, $8 * 100}' > ${tissue}_5mC.bedGraph &&
    
    echo "chrom start end pct_met cnt_C cnt_T cnt_N" | tr ' ' '\t' > 5mC_vs_alt_${tissue}.bed &&
    intersectBed -sorted -a ${tissue}_5mC.bedGraph -b ${tissue}_DNA.chr18.pct_alt.bedGraph -wa -wb \
    | groupBy -g 1,2,3,4 -c 8,9,10 -o sum,sum,sum >> 5mC_vs_alt_${tissue}.bed &&
     
    rm ${tissue}_5mC.bedGraph &&
    rm ${tissue}_DNA.chr18.pct_alt.bedGraph
done

tableCat.py -H -id tissue -i 5mC_vs_alt_*.bed -r '5mC_vs_alt_' | sed 's/.bed//' > 5mC_vs_alt.bed
rm 5mC_vs_alt_*.bed
```

### Plot 5mC vs alternative allele frequency

```
R
library(data.table)
library(scales)
library(ggplot2)

bdg<- fread('5mC_vs_alt.bed')
bdg[, false_pct_umet := 100 * cnt_T / (cnt_C + cnt_T)]
bdg[, pct_alt := 100 * (cnt_T + cnt_N) / (cnt_C + cnt_T + cnt_N)]
bdg[, pct_met_adj := pct_met + false_pct_umet]

bdg$cnt_tot<- rowSums(bdg[, list(cnt_C, cnt_T, cnt_N)])
gg<- ggplot(data= bdg[end < 40*1e6], aes(x= pct_alt)) +
    geom_histogram(binwidth= 1) +
    facet_wrap(~tissue)
    
metMarAltTum<- merge(
    bdg[tissue == 'margin', list(chrom, start, end, metMargin= pct_met_adj)],
    bdg[tissue == 'tumor', list(chrom, start, end, altTumor= pct_alt)],
by= c('chrom', 'start', 'end'))

gg<- ggplot(data= metMarAltTum[seq(1, nrow(metMarAltTum), length.out= 50000)], aes(x= metMargin, y= altTumor)) +
    geom_point(size= 0.5, colour= 'grey30', alpha= 0.20) +
    xlab('% 5mC margin') +
    ylab('% Alt. allele tumor') +
    ggtitle('Methylation and mutation frequency at CpG sites')
ggsave('5mCMargin_AltAllTumor.pdf', w= 12/2.54, h= 12/2.54)
system('rsync --remove-source-files 5mCMargin_AltAllTumor.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/')

gg<- ggplot(data= metMarAltTum[seq(1, nrow(metMarAltTum), length.out= 50000)], aes(x= metMargin, y= altTumor)) +
    geom_point(size= 0.5, colour= 'grey30', alpha= 0.20) +
    geom_smooth() + 
    xlab('% 5mC margin') +
    ylab('% Alt. allele tumor') +
    ylim(0, 10) +
    ggtitle('Methylation and mutation frequency at CpG sites')
ggsave('5mCMargin_AltAllTumor-zoom.pdf', w= 12/2.54, h= 12/2.54)
system('rsync --remove-source-files 5mCMargin_AltAllTumor-zoom.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/')

metMarAltTum[, metMarginTier := ifelse(metMargin < 10, 'low', ifelse(metMargin < 50, 'mid', 'high'))]
gg<- ggplot(data= metMarAltTum, aes(x= altTumor)) +
    geom_histogram(aes(y = ..density..)) +
    facet_wrap(~metMarginTier)
boxplot(altTumor ~ metMarginTier, data= metMarAltTum, outline= FALSE, notch= TRUE)

bdg[, pct_bin := ifelse(pct_met_adj < 10, '< 10%', '> 10%')]
mdn<- bdg[!is.na(pct_met_adj), list(median= median(pct_alt), avg= mean(pct_alt)), by= pct_bin]
gg<- ggplot(data= bdg[!is.na(pct_met_adj)], aes(x= pct_bin, y= ifelse(pct_alt > 10, 10, pct_alt))) +
    geom_violin(scale= 'area') +
    geom_point(data= mdn, aes(x= pct_bin, y= median), shape= 95, size= 5) + 
    xlab('5mC level') +
    ylab('% alternative allele') +
    scale_y_continuous(breaks=pretty_breaks(n=6)) +
    ggtitle('Relationship between methylation and mutation frequency') +
    facet_wrap(~tissue)
ggsave('5mC_vs_mut_freq.pdf', w= 20, h= 12, units= 'cm')
system('rsync --remove-source-files 5mC_vs_mut_freq.pdf $mac_office:~/Tritume/')

keep<- seq(1, nrow(bdg), length.out= 50000)
gg<- ggplot(data= bdg[keep], aes(x= pct_met_adj, y= pct_alt)) +
    geom_point(alpha= 0.2, size= 0.75) +
    xlim(0, 100) +
    ylim(0, 60) +
    facet_wrap(~tissue)
ggsave('5mC_vs_mut_freq-scatter.pdf', w= 24, h= 14, units= 'cm')
system('rsync --remove-source-files 5mC_vs_mut_freq-scatter.pdf $mac_office:~/Tritume/')

## Tumor-Margin difference in methylation and mutation frequency
pct_met_adj<- dcast.data.table(data= bdg, chrom + start ~ tissue, value.var= 'pct_met_adj')
pct_met_adj[, pct_met_diff := tumor - margin]

pct_alt<- dcast.data.table(data= bdg, chrom + start ~ tissue, value.var= 'pct_alt')
pct_alt[, pct_alt_diff := tumor - margin]
tmDiff<- merge(pct_met_adj, pct_alt, by= c('chrom', 'start'), suffixes= c('.met', '.alt'))
tmDiff<- tmDiff[complete.cases(tmDiff)]

keep<- seq(1, nrow(tmDiff), length.out= 120000)
gg<- ggplot(data= tmDiff[keep], aes(x= pct_met_diff, y= pct_alt_diff)) +
    geom_point(alpha= 0.2, size= 1) +
    xlab('5mC [Tumor - Margin]') +
    geom_vline(x= 0, colour= 'red', linetype= 'dashed', size= 0.5) +
    geom_hline(y= 0, colour= 'red', linetype= 'dashed', size= 0.5) + 
    ylab('Alt allele [Tumor - Margin]') +
    ggtitle('Relationship between methylationa and\nalternative allele frequency')
ggsave('5mC_vs_mut_tumor-margin.pdf', w= 24, h= 14, units= 'cm')
system('rsync --remove-source-files 5mC_vs_mut_tumor-margin.pdf $mac_office:~/Tritume/')
```

