## Segments and genomic features

We want to see the spatial relationship between HMM segments and genomic features (exons, introns, UTRs etc).

Possible approaches:

* `bedtools reldist`
* `bedtools jaccard`
* gat
* deepTools

#### Get features to analyze

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151124_5mC_TvsM_features/

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
grep 'chr18' /nas/sblab_data1/berald01/projects/20150501_methylation_brain/reference_data/$bed | cut -f1-3 > ${bed/hg19/chr18}.tmp
done
```

#### Segments to annotate

Collect states where M > T or vice versa, leave out state "2" from the 3 state HMM.

```
awk '$4 == 1' /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/oxBS_TvsM.states.bedGraph > oxBS_TvsM.T.bed
awk '$4 == -1' /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/oxBS_TvsM.states.bedGraph > oxBS_TvsM.M.bed
```

#### Test for association

Using [gat](http://gat.readthedocs.org/en/latest/):

```
## Workspace
grep 'chr18' /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.bed > chr18.size.bed

for x in chr18*.bed.tmp
do
    echo "gat-run.py -w chr18.size.bed -a $x -s oxBS_TvsM.T.bed -n 10000 > $x.gat.T.bed" > ${x}.T.tmp.sh
    echo "gat-run.py -w chr18.size.bed -a $x -s oxBS_TvsM.M.bed -n 10000 > $x.gat.M.bed" > ${x}.M.tmp.sh
done
ls *.tmp.sh | xargs --max-procs=24 --max-args=1 bash && rm *.tmp.sh
tableCat.py -i *.gat.*.bed -r '.*gat\.' | grep -v '#' > /tmp/gat.out.txt
(head -n 1 /tmp/gat.out.txt; grep -v 'track' /tmp/gat.out.txt) > gat.out.txt

rm /tmp/gat.out.txt
rm *.tmp.*
```

#### Plot association

```R
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151124_5mC_TvsM_features/
R
library(data.table)
library(scales)
library(ggplot2)

gat<- fread('gat.out.txt')
setnames(gat, 'M.bed', 'Tumor_is')
gat[, Tumor_is := ifelse(Tumor_is == 'M.bed', 'hypo_M', 'hyper_M')]
gat[, annotation := sub('.bed.tmp', '', annotation)]
gat[, annotation := sub('.*\\.', '', annotation)]
gat[, annotation := sub('promoters_1000upstream_0downstream', 'promoter', annotation)]

gat<- gat[order(Tumor_is, -l2fold)]
gat[, annotation := factor(annotation, levels= unique(annotation))]

gg<- ggplot(data= gat, aes(x= annotation, y= l2fold, by= Tumor_is, fill= Tumor_is)) +
    geom_bar(stat= 'identity', width = 0.6, position= position_dodge(width=0.6)) +
    theme(axis.text.x= element_text(angle= 45, hjust= 1)) +
    scale_y_continuous(breaks= pretty_breaks(10)) +
    ylab('Enrichment (log2FC)') +
    xlab('') +
    ggtitle('Association between genomic features and methylated blocks')
ggsave('barplot.gat.pdf', w= 20, h= 14, units= 'cm')
system('rsync --remove-source-files barplot.gat.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151124_5mC_TvsM_features/')

```


