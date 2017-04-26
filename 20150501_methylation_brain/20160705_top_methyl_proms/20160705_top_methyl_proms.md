<!-- MarkdownTOC -->

- [Top \(hydroxy\)methylated promoters and gene expression](#top-hydroxymethylated-promoters-and-gene-expression)

<!-- /MarkdownTOC -->

Top (hydroxy)methylated promoters and gene expression
=====================================================

We want to investigate some properties of the top most methylated and hydroxymethylated promoters.

*E.g.* 

* We need to assign to the promoters table the number of CpG and depth in order to ignore promoters with little evidence of methylation
* Pick the top 5hmC promoters or differential promoters between tumor and margin
* Assign gene expression
* Plot a bars or else...

Get data files of 5(h)mC in promoters and gene expression

```
cd /nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160705_top_methyl_proms/

scp uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/tsgenes_oncogenes/data/methylation_cpg_promoter/*.cpg_genes.promoters.1000.sorted.bed tmp/

R
library(data.table)
library(ggplot2)

## Promoter 5(h)mC
proms<- fread('tableCat.py -r "\\..*" -i tmp/ear*.cpg_genes.promoters.1000.sorted.bed')
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
               tum_mc= ear045_T3oxBS ,
               tum_hmc= ear044_T3BS - ear045_T3oxBS)]
prommet[, mar_hmc := ifelse(mar_hmc < 0, 0, mar_hmc)]
prommet[, tum_hmc := ifelse(tum_hmc < 0, 0, tum_hmc)]

## Table Gene Expression
tx<- fread('/mnt/nfs/nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq/tx_quant_lfc.bed')

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
    ggtitle('Tumor - Margin difference: Expression vs promoter 5mC')
ggsave('tum-mar_mc_vs_expr.png', width= 14, height= 12, units= 'cm')
system('rsync --remove-source-files tum-mar_mc_vs_expr.png $mac_office:~/git_sblab/projects/trunk/20150501_methylation_brain/20160705_top_methyl_proms/figures/')
```

<img src=figures/tum-mar_mc_vs_expr.png width="600">


* Tumor - Margin difference

```
gg<- ggplot(data= metexpr, aes(x= tum_hmc - mar_hmc, y= tpmLog2FC, 
        col= densCols(metexpr$tum_hmc - metexpr$mar_hmc, metexpr$tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))) + 
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    ggtitle('Tumor - Margin difference: Expression vs promoter 5hmC')
ggsave('tum-mar_hmc_vs_expr.png', width= 14, height= 12, units= 'cm')
system('rsync --remove-source-files tum-mar_hmc_vs_expr.png $mac_office:~/git_sblab/projects/trunk/20150501_methylation_brain/20160705_top_methyl_proms/figures/')
```
<img src=figures/tum-mar_hmc_vs_expr.png width=600>

* hmC relative conc
```
gg<- ggplot(data= metexpr, aes(x= (tum_hmc/tum_mc) - (mar_hmc/mar_mc), y= tpmLog2FC, 
    col= densCols((tum_hmc/tum_mc) - (mar_hmc/mar_mc), tpmLog2FC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))) + 
    geom_point(size= 0.1) +
    geom_hline(yintercept= 0, col= 'black', linetype= 'dashed') +
    geom_vline(xintercept= 0, col= 'black', linetype= 'dashed') +
    scale_color_identity() +
    ggtitle('Tumor - Margin difference: Expression vs rel.conc. 5hmC')
ggsave('tum-mar_rel_hmc_vs_expr.png', width= 14, height= 12, units= 'cm')
system('rsync --remove-source-files tum-mar_rel_hmc_vs_expr.png $mac_office:~/git_sblab/projects/trunk/20150501_methylation_brain/20160705_top_methyl_proms/figures/')
```
<img src=figures/tum-mar_rel_hmc_vs_expr.png width=600>

```
## Annotate promoters
library(biomaRt)
## Hosts can be found at http://www.ensembl.org/info/website/archives/index.html
mart<- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "dec2015.archive.ensembl.org", dataset= 'hsapiens_gene_ensembl')
genedescr<- data.table(getBM(
    filters= "external_gene_name",
    attributes= c("external_gene_name", "description"),
    values= metexpr$gene_name, 
    mart=mart))
genedescr<- genedescr[gsub("^//s+|//s+$", "", description) != ""] ## Remove genes w/o descr
genedescr<- genedescr[, list(description= paste(description, collapse= "|")), by= external_gene_name] ## Collapse dup descr
stopifnot(length(genedescr$external_gene_name) == length(unique(genedescr$external_gene_name)))

## Write-out table of promoter methyl, expr, annotation
xmrg<- merge(metexpr, genedescr, by.x= 'gene_name', by.y= 'external_gene_name', all.x= TRUE)
write.table(xmrg[, list(
    gene_name,
    mar_mc,
    mar_hmc,
    tum_mc,
    tum_hmc,
    tum_vs_mar_mc= tum_mc - mar_mc,
    tum_vs_mar_hmc= tum_hmc - mar_hmc,
    tum_vs_mar_rel_hmc= (tum_hmc/tum_mc) - (mar_hmc/mar_mc),
    tpmLog2FC,
    tpmLog2Avg,
    description
    )], 'genex_vs_methyl.txt', sep= '\t', row.names= FALSE, quote= FALSE)
system('gzip genex_vs_methyl.txt')
system('rsync --remove-source-files genex_vs_methyl.txt.gz $mac_office:~/git_sblab/projects/trunk/20150501_methylation_brain/20160705_top_methyl_proms/data/')
```

