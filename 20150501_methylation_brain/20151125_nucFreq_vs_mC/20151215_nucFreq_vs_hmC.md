## Nucleotide frequency vs methylation

Same as [20151125 nucFreq vs mC](20151125_nucFreq_vs_mC.md).

We want to see if there is a correlation between nucleotide frequencies and 5hmC at CpG sites.

For this we need:

* Table of 5hmC. See [20151102 Tumor vs Margin by bayesian approach](20151102_Tumor_vs_Margin_by_bayesian_approach)
* Table of nucleotide freqs from Sergio. See [20151124 SNV analysis](20151124_SNV_analysis)

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/

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

    ## 
    bigWigToBedGraph -chrom=chr18 ../20151106_genome_methyl/${tissue}.pct_5hmC.bw  /dev/stdout \
    | intersectBed -sorted -a - -b ${tissue}_DNA.chr18.pct_alt.bedGraph -wa -wb \
    | groupBy -g 1,2,3,4 -c 8,9,10 -o sum,sum,sum | less
    
    ## Hydroxymethylation
    bioawk -tc hdr 'NR > 1 {print $chrom, $start, $end, $mode * 100}' ../20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/hmC_${tissue}_posterior.txt.gz > ${tissue}_5hmC.bedGraph &&
    
    echo "chrom start end pct_met cnt_C cnt_T cnt_N" | tr ' ' '\t' > 5hmC_vs_alt_${tissue}.bed &&
    intersectBed -sorted -a ${tissue}_5hmC.bedGraph -b ${tissue}_DNA.chr18.pct_alt.bedGraph -wa -wb \
    | groupBy -g 1,2,3,4 -c 8,9,10 -o sum,sum,sum >> 5hmC_vs_alt_${tissue}.bed &&
     
    rm ${tissue}_5hmC.bedGraph &&
    rm ${tissue}_DNA.chr18.pct_alt.bedGraph
done

tableCat.py -H -id tissue -i 5hmC_vs_alt_*.bed -r '5hmC_vs_alt_' | sed 's/.bed//' > 5hmC_vs_alt.bed
rm 5hmC_vs_alt_*.bed
```

### Plot 5hmC vs alternative allele frequency

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

metMarAltTum<- merge(
    bdg[tissue == 'margin', list(chrom, start, end, altMargin= pct_alt, metMargin= pct_met_adj)],
    bdg[tissue == 'tumor', list(chrom, start, end, altTumor= pct_alt)],
by= c('chrom', 'start', 'end'))

gg<- ggplot(data= metMarAltTum[seq(1, nrow(metMarAltTum), length.out= 50000)], aes(x= metMargin, y= altTumor)) +
    geom_point(size= 0.5, colour= 'grey30', alpha= 0.20) +
    xlab('% 5hmC margin') +
    ylab('% Alt. allele tumor') +
    ggtitle('Hydroxy-methylation and mutation frequency at CpG sites')
ggsave('5hmCMargin_AltAllTumor.pdf', w= 12/2.54, h= 12/2.54)
system('rsync --remove-source-files 5hmCMargin_AltAllTumor.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/')

dat<- metMarAltTum[altTumor > 0 & altMargin < 1][seq(1, nrow(dat), length.out= 20000)]
gg<- ggplot(data= dat, aes(x= ifelse(metMargin < 0, 0, metMargin), y= altTumor)) +
    geom_point(size= 0.5, colour= 'grey30', alpha= 0.20) +
    geom_smooth() + 
    xlab('% 5hmC margin') +
    ylab('% Alt. allele tumor') +
    ylim(0, 10) +
    ggtitle('Methylation and mutation frequency at CpG sites')
ggsave('5hmCMargin_AltAllTumor-zoom.pdf', w= 12/2.54, h= 12/2.54)
system('rsync --remove-source-files 5hmCMargin_AltAllTumor-zoom.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151125_nucFreq_vs_mC/')

a<- nrow(metMarAltTum[altMargin < 1 & metMargin < 20 & altTumor < 1])
b<- nrow(metMarAltTum[altMargin < 1 & metMargin < 20 & altTumor >= 1])
c<- nrow(metMarAltTum[altMargin < 1 & metMargin >= 20 & altTumor < 1])
d<- nrow(metMarAltTum[altMargin < 1 & metMargin >= 20 & altTumor >= 1])

fisher.test(matrix(c(a,b,c,d), ncol= 2, byrow= TRUE))

b/(a+b)
d/(c+d)

#### NOT USED ==================================================================

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

