

## 5(h)mC in promoters

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151211_mC_target_genes/

awk -v OFS='\t' '$3 == "exon" {tx= gensub(/.*transcript_id \"/, "", $0); gsub(/\".*/, "", tx);
    gene= gensub(/.*gene_id \"/, "", $0); gsub(/\".*/, "", gene); print $1, $4, $5, $7, tx, gene}' \
    /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > genes.bed

## TSS of each transcipt on + extended to promoter
grep -P '\t\+\t' genes.bed | sort -k5,5 -k1,1n -k2,2n -k3,3n | groupBy -g 5 -c 1,2,6 -o distinct,min,distinct \
| sort -k4,4 -k3,3n | groupBy -g 3,4 -c 1,2 -o distinct,distinct | awk -v OFS='\t' '{print $4, $1-1, $1, $2, $3, "+"}' | sort -k1,1 -k2,2n > tss.plus.tmp.bed
## - strand
grep -P '\t\-\t' genes.bed | sort -k5,5 -k1,1n -k2,2n -k3,3n | groupBy -g 5 -c 1,3,6 -o distinct,max,distinct \
| sort -k4,4 -k3,3n | groupBy -g 3,4 -c 1,2 -o distinct,distinct | awk -v OFS='\t' '{print $4, $1-1, $1, $2, $3, "-"}' | sort -k1,1 -k2,2n > tss.minus.tmp.bed

sort -k1,1 -k2,2n -k3,3n tss.minus.tmp.bed tss.plus.tmp.bed > promoters.bed

for bw in /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/*.bw
do
    out=`basename $bw .bw`
    nohup computeMatrix reference-point -p 8 -S $bw -R promoters.bed -a 5000 -b 5000 -out /dev/null --outFileNameData ${out}.tss.tab &
done

tableCat.py -S 1 -i *.tss.tab -r '.tss.tab' | cut -f 1,2,3,5 > profile.tss.txt
```

Plot profile. Checked that forward and reverse features are correctly handled by computeMatrix. It probably requires strand info on
column 6, which is the case here.

```R
R
library(data.table)
library(ggplot2)
tab<- fread('profile.tss.txt')
setnames(tab, names(tab), c('bin', 'avg', 'stdev', 'library_id'))

tab[ , sample := sub('\\..*', '', library_id)]
tab[ , modification := sub('.*pct_', '', library_id)]

ggplot(data= tab, aes(x= bin, y= avg)) +
    geom_line(aes(linetype= modification, color= sample), size= 0.5) +
    xlab('Transcription start site') +
    ylab('% 5(h)mC') +
    ggtitle('Modification profile around TSS')
ggsave('mod_profileAllGenes_TSS.pdf', w= 16/2.54, h= 8/2.54)
system('rsync --remove-source-files mod_profileAllGenes_TSS.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151211_mC_target_genes/')
```

## 5(h)mC in target genes

We want to see if genes involved in DNA modification (_e.g._ DNMT, TET) have some methylation pattern in
tumor vs margin.

List of genes prepared by Gordon (see email 30/11/2015). Quote gene names to make it easier to extract them from
GTF file.

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151211_mC_target_genes/

grep -f geneListQ.txt /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > geneList.gtf

## Find genes with no match in GTF 
grep -o -f geneListQ.txt /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf | sort | uniq > /tmp/selected.txt
grep -f /tmp/selected.txt -v geneListQ.txt
    "AID"
    "MBDL3"
    "MBDL4"
    "MBDL5"
    "CRSP2"
rm /tmp/selected.txt

# Make as bed file
awk -v OFS='\t' '$3 == "exon" {tx= gensub(/.*transcript_id \"/, "", $0); gsub(/\".*/, "", tx); gene= gensub(/.*gene_id \"/, "", $0); gsub(/\".*/, "", gene); print $1, $4, $5, $7, tx, gene}' geneList.gtf > geneList.bed

## TSS of each transcipt on + extended to promoter
grep -P '\t\+\t' geneList.bed | sort -k5,5 -k1,1n -k2,2n -k3,3n | groupBy -g 5 -c 1,2,6 -o distinct,min,distinct \
| sort -k4,4 -k3,3n | groupBy -g 3,4 -c 1,2 -o distinct,distinct | awk -v OFS='\t' '{print $4, $1-1, $1, $2, $3, "+"}' | sort -k1,1 -k2,2n > tss.plus.tmp.bed
## - strand
grep -P '\t\-\t' geneList.bed | sort -k5,5 -k1,1n -k2,2n -k3,3n | groupBy -g 5 -c 1,3,6 -o distinct,max,distinct \
| sort -k4,4 -k3,3n | groupBy -g 3,4 -c 1,2 -o distinct,distinct | awk -v OFS='\t' '{print $4, $1-1, $1, $2, $3, "-"}' | sort -k1,1 -k2,2n > tss.minus.tmp.bed

sort -k1,1 -k2,2n -k3,3n tss.minus.tmp.bed tss.plus.tmp.bed > targetPromoters.bed

for bw in /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/*.bw
do
    out=`basename $bw .bw`
    nohup computeMatrix reference-point -p 8 -S $bw -R targetPromoters.bed -a 5000 -b 5000 -out ${out}.targets.tss.mat.gz --outFileNameData ${out}.targets.tss.tab &
done

tableCat.py -S 1 -i *.targets.tss.tab -r '.tss.tab' | cut -f 1,2,3,5 > profile.targets.tss.txt
```

Plot profile as above.

```R
R
library(data.table)
library(ggplot2)
tab<- fread('profile.targets.tss.txt')
setnames(tab, names(tab), c('bin', 'avg', 'stdev', 'library_id'))

tab[ , sample := sub('\\..*', '', library_id)]
tab[ , modification := sub('.*pct_', '', library_id)]

ggplot(data= tab, aes(x= bin, y= avg)) +
    geom_line(aes(linetype= modification, color= sample), size= 0.5) +
    xlab('Transcription start site') +
    ylab('% 5(h)mC') +
    ggtitle('Modification profile around TSS for genes involved in methylation metabolism')
ggsave('mod_profileTargetGenes_TSS.pdf', w= 16/2.54, h= 8/2.54)
system('rsync --remove-source-files mod_profileTargetGenes_TSS.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151211_mC_target_genes/')
```

### Plot individual genes

```R
R
library(data.table)
library(ggplot2)
mat<- read.table('tumor.pct_5mC.targets.tss.mat.gz', skip= 1)
names(mat)<- c('chrom', 'start', 'end', 'gene_name', 'score', 'strand', seq(-4990, 5000, by= 10))

pdf('profilesTargets.pdf')
for(i in 1:nrow(mat)){
    print(i)
    y<- as.numeric(mat[i, 7:ncol(mat)])
    naSkip<- !is.na(y)
    y<- y[naSkip]
    x<- seq(-4990, 5000, by= 10)[naSkip]
    gg<- ggplot(data= NULL, aes(x= x, y= y)) +
        geom_point() +
        geom_line()
    print(gg)
}
dev.off()
system('rsync --remove-source-files profilesTargets.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151211_mC_target_genes/')
```

chr1	32756708	32757708	HDAC1	NM_004964	+
chr1	44678125	44679125	DMAP1	NM_001034023,NM_001034024,NM_019100	+
chr1	68149860	68150860	GADD45A	NM_001199741,NM_001199742,NM_001924	+
chr10	69643427	69644427	SIRT1	NM_012238	+
chr10	69643939	69644939	SIRT1	NM_001142498	+
chr10	70319117	70320117	TET1	NM_030625	+
chr10	89622195	89623195	PTEN	NM_000314	+
chr10	131264454	131265454	MGMT	NM_002412	+
chr12	53772979	53773979	SP1	NM_001251825,NM_138473	+
chr12	53773428	53774428	SP1	NM_003109	+
chr12	57915659	57916659	MBD6	NM_052897	+
chr12	104358593	104359593	TDG	NM_003211	+
chr12	109534399	109535399	UNG	NM_080911	+
chr12	109534923	109535923	UNG	NM_003362	+
chr13	32888617	32889617	BRCA2	NM_000059	+
chr15	78440719	78441719	IDH3A	NM_005530	+
chr16	67595310	67596310	CTCF	NM_001191022,NM_006565	+
chr17	80476594	80477594	FOXK2	NM_004514	+
chr19	7048351	7049351	MBD3L2	NM_144614	+
chr2	74272450	74273450	TET3	NM_144993	+
chr2	148777580	148778580	MBD5	NM_018328	+
chr20	31349191	31350191	DNMT3B	NM_001207055,NM_001207056,NM_006892,NM_175848,NM_175849	+
chr20	31366658	31367658	DNMT3B	NM_175850	+
chr22	39352527	39353527	APOBEC3A	NM_001193289,NM_145699	+
chr22	39377405	39378405	APOBEC3B	NM_004900	+
chr22	39416118	39417118	APOBEC3D	NM_152426	+
chr22	39435673	39436673	APOBEC3F	NM_001006666,NM_145298	+
chr22	39492229	39493229	APOBEC3H	NM_001166002,NM_001166003,NM_001166004,NM_181773	+
chr22	41600313	41601313	L3MBTL2	NM_031488	+
chr3	178865311	178866311	PIK3CA	NM_006218	+
chr4	106066032	106067032	TET2	NM_017628	+
chr4	106066842	106067842	TET2	NM_001127208	+
chr6	28316691	28317691	ZKSCAN3	NM_001242894,NM_001242895,NM_024493	+
chr6	41019940	41020940	APOBEC2	NM_006789	+
chr6	41513164	41514164	FOXP4	NM_001012426,NM_001012427,NM_138457	+
chr7	18125572	18126572	HDAC9	NM_001204144	+
chr7	18534369	18535369	HDAC9	NM_001204145,NM_001204146,NM_014707,NM_178423	+
chr7	18534885	18535885	HDAC9	NM_058176,NM_178425	+
chr7	18547900	18548900	HDAC9	NM_001204147,NM_001204148	+
chr7	99646417	99647417	ZSCAN21	NM_145914	+
chr9	4984245	4985245	JAK2	NM_004972	+
chr9	131444934	131445934	SET	NM_001122821	+
chr9	131446365	131447365	SET	NM_001248000	+
chr9	131450509	131451509	SET	NM_003011	+
chr9	131451254	131452254	SET	NM_001248001	+
chrX	48659487	48660487	HDAC6	NM_006044	+
chrX	70751912	70752912	OGT	NM_181672,NM_181673	+
