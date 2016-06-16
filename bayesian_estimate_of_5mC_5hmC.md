<!-- MarkdownTOC -->

- [Estimating tumor vs margin differences via Bayesian approach](#estimating-tumor-vs-margin-differences-via-bayesian-approach)
    - [Data preparation](#data-preparation)
    - [Margin, high and low](#margin-high-and-low)
    - [Tumor, high and low](#tumor-high-and-low)
    - [Apply simPostDiff to all CpGs:](#apply-simpostdiff-to-all-cpgs)
    - [Convert to data.table](#convert-to-datatable)
    - [You might want to quit R and reload this obj to clean up the memory](#you-might-want-to-quit-r-and-reload-this-obj-to-clean-up-the-memory)
    - [Add to posterior datatable the position and raw counts](#add-to-posterior-datatable-the-position-and-raw-counts)
    - [Write out tables](#write-out-tables)
- [Estimating methylation and hydroxymethylation at base resolution](#estimating-methylation-and-hydroxymethylation-at-base-resolution)
- [Bayesian estimate of 5mC in tumor and margin](#bayesian-estimate-of-5mc-in-tumor-and-margin)
- [Bayesian estimate of differential methylation](#bayesian-estimate-of-differential-methylation)

<!-- /MarkdownTOC -->

Estimating tumor vs margin differences via Bayesian approach
============================================================

<!-- 
Working dir and input files

Methylation files taken from Sergio. see `20150501_methylation_brain/20160209_regenerating_methylation_files.md`

```
cd /nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160613_bayes/
scp uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/ear04*BS.cpg.bedGraph.gz ./

-->

Collapse methylation and total count at CpG sites
-------------------------------------------------

In order to increase statistical power, the counts of methylated and total reads at 
C and G at each CpG site were collapse. In other words, each CpG site will be treated as a single unit
rather than a C and G pair. 

```
intersectBed -a hg19.allCpG.bed.gz -b $bdg -wa -wb -sorted \
| awk -v OFS="\t" '{print $1, $2, $3, $12, $13}' \
| groupBy -i - -g 1,2,3 -c 4,5 -o sum,sum \
| awk -v OFS="\t" '{if($7 > 1) print $2, $3, $4, $5, $6}' > $out

intersectBed -a hg19.allCpG.bed.gz -b $bdg -wa -wb -sorted \
| awk -v OFS="\t" '{print $1, $2, $3, $12, $13, $4}' \
| groupBy -i - -g 6 -c 1,2,3,4,5,6 -o first,min,max,sum,sum,count \
| awk -v OFS="\t" '{if($7 > 1) print $2, $3, $4, $5, $6}' > $out
```

File `hg19.allCpG.bed.gz` is a bed file of CpG positions in the hg19 genome. It must contain 7 columns although only the first four 
are relevant. This file can be produced with the script [fastaRegexFinder.py](https://github.com/dariober/bioinformatics-cafe/blob/master/fastaRegexFinder.py)
as follows:

```
fastaRegexFinder.py -f  genome.fa -r CG --noreverse | gzip > hg19.allCpG.bed.gz
```

`$bdg` is the bed file output of `bam2methylation.py`. 
It must contain in the 5th and 6th column the methylated and total read counts, respectively.


Difference in methylation between tumor and margin
--------------------------------------------------

This part aims at estimating differences in methylation between the tumor and margin sample. 
Therefore, only the oxBS data is considered.

```
R
library(ggplot2)
library(data.table)
library(gridExtra)
library(parallel)
library(fitdistrplus)

Mode <- function(x) {
  ## Get mode of vector x. I.e. find the max of the kernel density of x
  ## Default bandwidth is not the best, but it's fast
  z<- density(x)    
  return( z$x[z$y==max(z$y)] )
}

## Data preparation
tum<- fread('zcat ear045_T3oxBS.cpg.bedGraph.gz')
mar<- fread('zcat ear043_M8oxBS.cpg.bedGraph.gz')
stopifnot(all(mar$V3 - mar$V2 == 2))
stopifnot(all(tum$V3 - tum$V2 == 2))

mar[, V3 := NULL]
tum[, V3 := NULL]

xn<- c('chrom', 'start', 'cnt_met', 'cnt_tot')
setnames(mar, names(mar), xn)
setnames(tum, names(tum), xn)
bdg<- merge(mar, tum, by= c('chrom', 'start'), suffixes= c('_mar', '_tum'))
rm(mar, tum)
```

Methylation tends to follow a bimodal distribution: 

```
set.seed(1234)
xn<- sample(1:nrow(bdg), size= nrow(bdg)/100, replace= FALSE)
gg<- ggplot(data= bdg[xn]) +
    geom_histogram(aes(x= 100 * cnt_met_mar/cnt_tot_mar), fill= 'blue', alpha= 0.4) +
    geom_histogram(aes(x= 100 * cnt_met_tum/cnt_tot_tum), fill= 'red', alpha= 0.4) +
    ylab('N. CpG x1000') +
    xlab('% 5mC') +
    annotate("text", x = c(0, 0), y = c(40000, 38000), label = c("Margin", 'Tumor'), 
        colour= c('blue', 'red'), vjust= 1) +
    ggtitle('Distribution of 5mC in margin and tumor')
ggsave('hist_5mc_margin_tumor.png', w= 14, h= 12, units= 'cm')
```

<img src="figures/hist_5mc_margin_tumor.png" width="600">


### Choosing prior distribution for Bayesian estimation

In order to capture the bimodal distribution of 5mC, the prior distributions 
are empirically constructed from a mixture of two beta distributions. 
The `CUTOFF` parameter determines the boundary of the two distributions. This parameter is
chosen by visually inspecting the distribution of 5mC in tumor and margin.

`priorBetaParam_*_High` contains the two parameters of a beta distribution fitted to 
the datapoints > CUTOFF. Conversely `priorBetaParam_*_Low` contains the beta parameters
for the datapoints < CUTOFF.

```
CUTOFF<- 0.5

## Margin, high and low
pct_met<- bdg[cnt_met_mar/cnt_tot_mar > CUTOFF, cnt_met_mar/cnt_tot_mar]
priorBetaParam_M_High<- fitdist(pct_met[sample(1:length(pct_met), size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))
pct_met<- bdg[cnt_met_mar/cnt_tot_mar <= CUTOFF, cnt_met_mar/cnt_tot_mar]
priorBetaParam_M_Low<- fitdist(pct_met[sample(1:length(pct_met), size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

## Tumor, high and low
pct_met<- bdg[cnt_met_tum/cnt_tot_tum > CUTOFF, cnt_met_tum/cnt_tot_tum]
priorBetaParam_T_High<- fitdist(pct_met[sample(1:length(pct_met), size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))
pct_met<- bdg[cnt_met_tum/cnt_tot_tum <= CUTOFF, cnt_met_tum/cnt_tot_tum]
priorBetaParam_T_Low<- fitdist(pct_met[sample(1:length(pct_met), size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))
```

In order to plot and make use of the prior distributions we sample a number of data points from each pair (high and low).
The number of data points sampled from the high and low distributions is proportional to the number
of CpGs falling below and above the cutoff.

```
prop_M_high<- nrow(bdg[cnt_met_mar/cnt_tot_mar > CUTOFF]) / nrow(bdg)
prop_T_high<- nrow(bdg[cnt_met_tum/cnt_tot_tum > CUTOFF]) / nrow(bdg)

NSIM<- 100000
rndPrior_M<- c(
    rbeta(10 * NSIM * prop_M_high, priorBetaParam_M_High$estimate[1], priorBetaParam_M_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_M_high), priorBetaParam_M_Low$estimate[1], priorBetaParam_M_Low$estimate[2]))
rndPrior_T<- c(
    rbeta(10 * NSIM * prop_T_high, priorBetaParam_T_High$estimate[1], priorBetaParam_T_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_T_high), priorBetaParam_T_Low$estimate[1], priorBetaParam_T_Low$estimate[2]))

gg<- ggplot(data= NULL) +
    geom_histogram(aes(x= 100 * rndPrior_M, y= ..density..), fill= 'blue', alpha= 0.4) +
    geom_histogram(aes(x= 100 * rndPrior_T, y= ..density..), fill= 'red', alpha= 0.4) +
    xlab('% 5mC') + 
    ggtitle('Prior distributions of 5mC in margin and tumor')
ggsave('bimodalBetaPrior_5mc_margin_tumor.png', w= 14, h= 12, units= 'cm')
```

<img src="figures/bimodalBetaPrior_5mc_margin_tumor.png" width="600">

### Estimating the posterior distributions

For each CpG update the prior using the observed data.

```
simPostDiff<- function(x, priorBetaParam_T_High, priorBetaParam_T_Low, 
                          priorBetaParam_M_High, priorBetaParam_M_Low,
                          prop_T_high, prop_M_high, nsim= 20000){
    ## Update priors given data in x
    ## x: 
    ##      Vectors of length 4 giving the observed data: 
    ##      1) count methylated in margin 2) count total in margin
    ##      3) count methylated tumor 4) count total tumor
    ## priorBetaParam_T/M_High/Low:
    ##      Parameters of the beta distribution of the low and high mC tails, 
    ##      for tumor and margin
    ## prop_T/M_high: 
    ##      Proportion of the entire mC distribution belonging to the 'high' tail
    ## nsim:
    ##      Number of random samples to draw for simulation
    ## Returns:
    ##      List of 3 vectors: Quantiles and mode of the posteriors of: 
    ##      1) Tumor 5mC, 2) Margin 5mC, 3) difference T-M
    
    p<- c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995)
    
    x<- unlist(x)
    stopifnot(length(x) == 4)
    stopifnot(all(!is.na(x)))

    yM = x[1]   # Met in Margin
    nM = x[2]   # Counts
    yT = x[3]   # Met in tumor
    nT = x[4]

    ## Some simple sanity check
    stopifnot(yT <= nT)
    stopifnot(yM <= nM)

    post_T<- c(
        rbeta(nsim * prop_T_high, yT + priorBetaParam_T_High$estimate[1], (nT - yT) + priorBetaParam_T_High$estimate[2]),
        rbeta(nsim * (1-prop_T_high), yT + priorBetaParam_T_Low$estimate[1], (nT - yT) + priorBetaParam_T_Low$estimate[2]))
    post_M<- c(
        rbeta(nsim * prop_M_high, yM + priorBetaParam_M_High$estimate[1], (nM - yM) + priorBetaParam_M_High$estimate[2]),
        rbeta(nsim * (1-prop_M_high), yM + priorBetaParam_M_Low$estimate[1], (nM - yM) + priorBetaParam_M_Low$estimate[2]))

    nn<- sprintf('p%s', p)
    Tq<- quantile(post_T, p);
    names(Tq)<- nn
    Tq<- c(Tq, mode= Mode(post_T))
    Mq<- quantile(post_M, p)
    names(Mq)<- nn
    Mq<- c(Mq, mode= Mode(post_M))
    Dq<- quantile(post_T - post_M, p)
    names(Dq)<- nn
    Dq<- c(Dq, mode= Mode(post_T - post_M))
    
    qq<- list(post_T= Tq, post_M= Mq, post_D= Dq)
    return(qq)
}

## Apply simPostDiff to all CpGs:
datOxBS<- bdg[, list(cnt_met_mar,
                     cnt_tot_mar,
                     cnt_met_tum,
                     cnt_tot_tum)] ## Order of columns matters!
clus<- makeCluster(24)
clusterExport(clus, list('simPostDiff', 'Mode', 'priorBetaParam_T_High', 'priorBetaParam_T_Low', 'priorBetaParam_M_High', 'priorBetaParam_M_Low',
    'prop_T_high', 'prop_M_high'))
system.time({
outqq<- parRapply(clus, datOxBS, 
    function(x) simPostDiff(x, 
        priorBetaParam_T_High, priorBetaParam_T_Low, priorBetaParam_M_High, priorBetaParam_M_Low,
        prop_T_high, prop_M_high, nsim= 20000))
}) # MEMO: For NSIM=20000 elapsed time is ~53000s
stopCluster(clus)
save(outqq, file= 'outqq.tmp.Rdata')

## Convert to data.table
cnames<- names(outqq[[1]][[1]])
slots<- names(outqq[[1]])
postDT<- data.table(matrix(unlist(outqq), ncol= length(cnames) * length(slots), byrow= TRUE))
rm(outqq); gc()

xnames<-  apply(expand.grid(cnames, slots), 1, paste, collapse= '_')
setnames(postDT, names(postDT), xnames)
write.table(x= postDT, file= 'postDT.tmp.txt', row.names= FALSE, quote= FALSE)

## You might want to quit R and reload this obj to clean up the memory
postDT<- fread('postDT.tmp.txt')

## Add to posterior datatable the position and raw counts
postDT<- cbind(bdg, postDT)
rm(bdg); gc()

## Write out tables
write.table(x= 
  postDT[, list(chrom,
                start,
                cnt_met_mar,
                cnt_tot_mar,
                cnt_met_tum,
                cnt_tot_tum,
                p0.005= 100 * round(p0.005_post_D, 4),
                p0.025= 100 * round(p0.025_post_D, 4),
                p0.25=  100 * round(p0.25_post_D, 4),
                p0.5=   100 * round(p0.5_post_D, 4),
                p0.75=  100 * round(p0.75_post_D, 4),
                p0.975= 100 * round(p0.975_post_D, 4),
                p0.995= 100 * round(p0.995_post_D, 4),
                mode= 100 * round(mode_post_D, 4))],
  file= 'posterior_5mC_tum-mar.txt', row.names= FALSE, quote= FALSE, sep= '\t')

write.table(x= 
  postDT[, list(chrom,
                start,
                p0.005= 100 * round(p0.005_post_T, 4),
                p0.025= 100 * round(p0.025_post_T, 4),
                p0.25=  100 * round(p0.25_post_T, 4),
                p0.5=   100 * round(p0.5_post_T, 4),
                p0.75=  100 * round(p0.75_post_T, 4),
                p0.975= 100 * round(p0.975_post_T, 4),
                p0.995= 100 * round(p0.995_post_T, 4),
                mode= 100 * round(mode_post_T, 4))],
  file= 'posterior_5mC_tumor.txt', row.names= FALSE, quote= FALSE, sep= '\t')

write.table(x= 
  postDT[, list(chrom,
                start,
                p0.005= 100 * round(p0.005_post_M, 4),
                p0.025= 100 * round(p0.025_post_M, 4),
                p0.25=  100 * round(p0.25_post_M, 4),
                p0.5=   100 * round(p0.5_post_M, 4),
                p0.75=  100 * round(p0.75_post_M, 4),
                p0.975= 100 * round(p0.975_post_M, 4),
                p0.995= 100 * round(p0.995_post_M, 4),
                mode= 100 * round(mode_post_M, 4))],
  file= 'posterior_5mC_margin.txt', row.names= FALSE, quote= FALSE, sep= '\t')
system('gzip posterior_5mC*.txt')
```

How does the Bayesian estimate of Tumor-Margin compare to the simple difference between percentages (i.e. %tum - $mar)?

```
R
library(data.table)
library(scales)
library(ggplot2)

postDiff<- fread('zcat posterior_5mC_tum-mar.txt.gz')
xn<- seq(1, nrow(postDiff), length.out= 50000)
gg<- ggplot(data= postDiff[xn, ], aes(x= 100 * ((cnt_met_tum/cnt_tot_tum) - (cnt_met_mar/cnt_tot_mar)), y= mode,
    colour= ifelse(cnt_tot_mar + cnt_tot_tum > 200, 200, cnt_tot_mar + cnt_tot_tum))) +
    scale_colour_gradient2('Read count', low=muted("blue"), high=muted("red"), midpoint= 75, mid= muted('blue')) +
    geom_point(alpha= 0.50, size= 0.15) +
    geom_abline(intercept= 0, slope= 1, linetype= 'dashed') +
    xlab('Observed tumor-margin difference') +
    ylab('Bayesian posterior difference') +
    ggtitle('Tumor-Margin difference in 5mC\nObserved from raw counts vs Bayesian estimate')
ggsave('posterior_difference_5mc_margin_tumor.png', w= 14, h= 12, units= 'cm')

dat<- rbindlist(list(
    postDiff[, list(diff= mode, method= 'Posterior est.')],
    postDiff[, list(diff= 100 * ((cnt_met_tum/cnt_tot_tum) - (cnt_met_mar/cnt_tot_mar)), method= 'Obs. diff.')]
))
gg<- ggplot(data= dat, aes(x= diff, fill= method)) +
    geom_histogram(alpha= 0.6, colour= 'white', position= 'identity') +
    geom_vline(xintercept= 0, linetype= 'dashed') +
    xlab('5mC Tumor - Margin')
ggsave('posterior_difference_5mc_margin_tumor_hist.png', w= 14, h= 10, units= 'cm')
```
<img src="figures/posterior_difference_5mc_margin_tumor.png" width="500"> 
<img src="figures/posterior_difference_5mc_margin_tumor_hist.png" width="600">


Here we look at the size of the credibility intervals in relation to total read depth, i.e. depth in margin plus tumor.

```
postDiff[, depth_bin := round((cnt_tot_mar + cnt_tot_tum)/20)*20]
postDiff[, depth_bin := ifelse(depth_bin > 300, 300, depth_bin)]
postDiff[, depth_bin := factor(depth_bin, levels= sort(unique(depth_bin)))]
xn<- seq(1, nrow(postDiff), length.out= 500000)
gg<- ggplot(data= postDiff[xn,], aes(x= depth_bin, y= abs(p0.025 - p0.975))) +
    geom_boxplot(varwidth= TRUE, size= 0.2, outlier.size= 0.2) +
    xlab('Binned read count (margin + tumor)') +
    ylab('Width of CI [0.025 0.975]') +
    ggtitle('Width of credibility intervals in relation to read depth')
ggsave('credint_5mc_margin_tumor.png', w= 16, h= 10, units= 'cm')
```

<img src="figures/credint_5mc_margin_tumor.png" width="700"> 

<!-- 
PRO MEMORIA TO BE DELETED:

Estimating methylation and hydroxymethylation at base resolution
================================================================

It would be good to submit as supplementary material tables of:

* 5mC in tumor and margin, from oxBS. Columns would be CpG position, count C, count C+T, estimate of 5mC and significance of difference from 0 with baysian approach.
* As above but for 5hmC. There would be additional columns for counts in BS and oxBS. 
* Difference in 5mC between tumor and margin. Format would be the same as for 5hmC. Instead of BS vs oxBS use *oxBS_tumor* vs *oxBS_margin*
* Ideally: Difference in **5hmC** between tumor and margin. This would come from difference of differences (need to figure out how to do it in Bayesian framework)

Get most of this from https://github.com/sblab-bioinformatics/projects/tree/master/20150501_methylation_brain/20151102_Tumor_vs_Margin_by_bayesian_approach


Bayesian estimate of 5mC in tumor and margin
============================================



Bayesian estimate of differential methylation
=============================================
-->
