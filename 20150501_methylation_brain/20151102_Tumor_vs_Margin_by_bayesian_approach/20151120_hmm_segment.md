# Segmenting 5mC difference between tumor and margin

Use the raw difference, not the posterior from bayes approach. This is easier to handle and explain.

Prepare data by taking the difference of percentages and smoothing by rolling mean.

```R
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20160128_hmmSegment

R
library(RHmm)
library(zoo)
library(ggplot2)
margin<- fread('bigWigToBedGraph ../20151106_genome_methyl/margin.pct_5mC.bw /dev/stdout')
tumor<- fread('bigWigToBedGraph ../20151106_genome_methyl/tumor.pct_5mC.bw /dev/stdout')
tm<- merge(margin, tumor, by= c('V1', 'V2', 'V3'))
setnames(tm, names(tm), c('chrom', 'start', 'end', 'pct_mrg', 'pct_tum'))
tm[, diff := pct_mrg - pct_tum]
rmDiff<- tm[, list(rmDiff= rollmean(.SD$diff, 75, na.pad= TRUE)), by= chrom]
tm[, rmDiff := rmDiff$rmDiff]

tm<- tm[!is.na(rmDiff)]
write.table(tm[, list(chrom, start, end, rmDiff)], '/tmp/rmDiff.bedGraph', col.names= FALSE, row.names= FALSE, quote= FALSE, sep= '\t')

## Test HMM params
hmmfit<- HMMFit(tm$rmDiff[1:100000], nStates= 3, dis= 'NORMAL')
quit()
# ================
```

Segment difference

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20160128_hmmSegment

chroms=`cut -f 1 /tmp/rmDiff.bedGraph | sort | uniq`
for chrom in $chroms
do
echo "grep -w $chrom /tmp/rmDiff.bedGraph | segmentPvalueByHMM.R -i - -p 4 -n 3 > ${chrom}.tmp.hmm.txt" > ${chrom}.tmp.hmm.sh
done
ls *.tmp.hmm.sh | xargs -n 1 -P 12 bash && rm *.tmp.hmm.sh
tableCat.py -i chr*.tmp.hmm.txt -H | rev | cut -f 2- | rev > oxBS_TvsM.hmm.bed

## Sanity check. Chrom names are equal. Must print no lines:
paste /tmp/rmDiff.bedGraph <(tail -n+2 oxBS_TvsM.hmm.bed) | awk '$1 != $5'

## Add position of cpg
hdr=`head -n 1 oxBS_TvsM.hmm.bed`
echo "chrom start end yData state probState_1 probState_2 probState_3 state_id" | tr ' ' '\t' > /tmp/oxBS_TvsM.hmm.bed
paste /tmp/rmDiff.bedGraph <(tail -n+2 oxBS_TvsM.hmm.bed) | cut -f 1,2,3,4,7- >> /tmp/oxBS_TvsM.hmm.bed
pigz /tmp/oxBS_TvsM.hmm.bed
mv  /tmp/oxBS_TvsM.hmm.bed ./

## bed file of blocks
zcat oxBS_TvsM.hmm.bed.gz | tail -n+2 \
| groupBy -g 1,5 -c 2,3 -o min,max \
| awk -v OFS="\t" '{print $1, $3, $4, $2}' > oxBS_TvsM.states.bedGraph

rm chr*tmp.hmm.txt
rm /tmp/rmDiff.bedGraph

pigz oxBS_TvsM.hmm.bed
```

## Summarize segmentation

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20160128_hmmSegment/
R
library(data.table)
library(ggplot2)
hdr<- fread('zcat oxBS_TvsM.hmm.bed.gz | head -n 10')
hmm<- fread('zcat oxBS_TvsM.hmm.bed.gz | grep -w chr18')
setnames(hmm, names(hmm), names(hdr))

hmm$state<- as.factor(hmm$state)
xrle<- rle(as.vector(hmm$state))
id<- rep(1:length(xrle$values), times= xrle$length)
hmm[, state_id := paste(state, id, sep= '_')]
segs<- hmm[, list(start= min(start), end= max(end), TvsM= mean(yData), .N), by= list(state, state_id)]
segs[, size := end-start]
gg<- ggplot(data= segs, aes(x= state, y= ifelse(N > 500, 500, N), fill= state, colour= state, alpha= 0.4)) +
    geom_boxplot(outlier.size= 0.5, width= 0.5) +
    xlab('State: Tumor vs Margin') +
    ylab('N. CpG sites') +
    ggtitle('N. of CpG sites in segments different\nbetween tumor and margin') +
    theme(legend.position= 'none')
ggsave('boxplot_N_CpG_oxBS_TvsM.hmm.pdf', w= 12, h= 12, units= 'cm')
system('rsync --remove-source-files boxplot_N_CpG_oxBS_TvsM.hmm.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20160128_hmmSegment/')

options(scipen= 10)
gg<- ggplot(data= segs, aes(x= state, y= size, fill= state, colour= state, alpha= 0.3)) +
    geom_violin() +
    scale_y_log10() + 
    xlab('State: Tumor vs Margin') +
    ylab('Segment size (bp)') +
    ggtitle('Segment size from HMM\nin tumor vs margin') +
    theme(legend.position= 'none')
ggsave('violin_Size_oxBS_TvsM.hmm.pdf', w= 12, h= 12, units= 'cm')
system('rsync --remove-source-files violin_Size_oxBS_TvsM.hmm.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20160128_hmmSegment/')

## NB: Make difference as Tumor - Margin
gg<- ggplot(data= hmm, aes(x= state, y= -(yData), fill= state, colour= state, alpha= 0.3)) +
    geom_boxplot(outlier.size= 0.5, width= 0.5) +
    xlab('State: Tumor vs Margin') +
    ylab('5mC % Tumor - % Margin') +
    ggtitle('Difference in 5mC between states') +
    theme(legend.position= 'none') +
    geom_hline(yintercept= 0, linetype= 'dashed', colour= 'grey20')
ggsave('boxplot_5mC_oxBS_TvsM.hmm.pdf', w= 12, h= 12, units= 'cm')
system('rsync --remove-source-files boxplot_5mC_oxBS_TvsM.hmm.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20160128_hmmSegment/')
```
