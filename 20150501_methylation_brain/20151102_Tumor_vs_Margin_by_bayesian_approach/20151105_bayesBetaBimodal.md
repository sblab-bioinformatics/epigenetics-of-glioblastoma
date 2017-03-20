<!-- MarkdownTOC -->

- [13/06/2016 Estimating Tumor vs Margin differences via Bayesian approach](#13062016-estimating-tumor-vs-margin-differences-via-bayesian-approach)
    - [Difference in methylation between tumor and margin](#difference-in-methylation-between-tumor-and-margin)
        - [Choosing prior distribution for Bayesian estimation](#choosing-prior-distribution-for-bayesian-estimation)
    - [Let's see the difference with larger](#lets-see-the-difference-with-larger)
- [Estimating Tumor vs Margin differences via Bayesian approach](#estimating-tumor-vs-margin-differences-via-bayesian-approach)
    - [Summary plots](#summary-plots)

<!-- /MarkdownTOC -->

13/06/2016 Estimating Tumor vs Margin differences via Bayesian approach
=======================================================================

Same as below, on whole genome, and tidy-up. Ideally this will go as supp mat online

Methylation files taken from Sergio. see `20150501_methylation_brain/20160209_regenerating_methylation_files.md`

```
cd /nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160613_bayes/
scp uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/ear04*BS.cpg.bedGraph.gz ./
```

Difference in methylation between tumor and margin
--------------------------------------------------

This part aims at estimating differences in methylation between the tumor and margin sample. 
Therefore, only the oxBS data is considered.

```
cd /nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160613_bayes/

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

### Choosing prior distribution for Bayesian estimation

In order to capture the bimodal distribution of 5mC, the prior distributions 
are empirically constructed from a mixture of two beta distribution. 
The `cutoff` parameter determines the boundary of the two distributions. This parameter is
chosen by visually inspecting the distribution of 5mC in tumor and margin.

`priorBetaParam_*_High` contains the two parameters of a beta distribution fitted to 
the datapoints > cutoff. Conversely `priorBetaParam_*_Low` contains the beta parameters
for the datapoints < cutoff.

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
}) # MEMO: For NSIM=20000 elapsed time is 53226
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

<!-- 
OLD ANALYSIS 
-->


## Let's see the difference with larger
clus<- makeCluster(24)
clusterExport(clus, list('simPostDiff', 'Mode', 'priorBetaParam_T_High', 'priorBetaParam_T_Low', 'priorBetaParam_M_High', 'priorBetaParam_M_Low',
    'prop_T_high', 'prop_M_high'))
outqq2<- parRapply(clus, datOxBS[1:10000,], 
    function(x) simPostDiff(x, 
        priorBetaParam_T_High, priorBetaParam_T_Low, priorBetaParam_M_High, priorBetaParam_M_Low,
        prop_T_high, prop_M_high, nsim= 200000))
stopCluster(clus)

postDT<- data.table(do.call('rbind', lapply(outqq[1:10000], function(x) c(x[['post_M']], library_id= 'ear045_T3oxBS'))))
postDT2<- data.table(do.call('rbind', lapply(outqq2, function(x) c(x[['post_M']], library_id= 'ear045_T3oxBS'))))
plot(postDT2$mode, postDT$mode)
which((as.numeric(postDT$mode) - as.numeric(postDT2$mode)) > 0.05)

plot(datOxBS[1:10000, cnt_met_mar/cnt_tot_mar], postDT$mode)


# Estimating Tumor vs Margin differences via Bayesian approach

```R
R
setwd('/nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18')

library(ggplot2)
library(data.table)
library(gridExtra)
library(parallel)

Mode <- function(x) {
  x<- round(x, 3)
  ux <- unique(round(x, 3))
  ux[which.max(tabulate(match(x, ux)))]
}

inv.logit<- function(x){
    exp(x)/(1+exp(x))
}

NSIM<- 100000

## Data preparation
bdg<- fread('zcat hiseq201509.chr18.cpg.bdg.gz')
bdg[, pct_met := cnt_met/cnt_tot]

smryAll<- bdg[, list(avg= mean(pct_met, na.rm= TRUE),  median= median(pct_met, na.rm= TRUE), stdev= sd(pct_met, na.rm= TRUE)), by= library_id]
#> smryAll
#      library_id       avg    median     stdev
#1:   ear042_M8BS 0.7500937 0.8461538 0.2588769
#2: ear043_M8oxBS 0.5765433 0.6363636 0.2447375
#3:   ear044_T3BS 0.6983177 0.8095238 0.2737420
#4: ear045_T3oxBS 0.6850772 0.7959184 0.2731436

gg<- ggplot(data= bdg[library_id %in% list('ear043_M8oxBS', 'ear045_T3oxBS'),], aes(x= pct_met)) +
    geom_histogram(fill= 'grey60', colour= 'grey20') +
    facet_wrap(~library_id)
ggsave('hist_pct_met_global.pdf', w= 20, h= 12, units= 'cm')
system('rsync -v --remove-source-files hist_pct_met_global.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151030_betaBinEstimate/')
# Only >x
smry<- bdg[pct_met > 0.5, list(avg= mean(pct_met, na.rm= TRUE), stdev= sd(pct_met, na.rm= TRUE)), by= library_id]

# Cross tabulate the comparison of interest
cnt_met<- dcast.data.table(data= bdg[library_id %in% c('ear043_M8oxBS', 'ear045_T3oxBS')], chrom+start+end ~ library_id, value.var= 'cnt_met')
cnt_tot<- dcast.data.table(data= bdg[library_id %in% c('ear043_M8oxBS', 'ear045_T3oxBS')], chrom+start+end ~ library_id, value.var= 'cnt_tot')
cmp<- merge(cnt_met, cnt_tot, by= c('chrom', 'start', 'end'), suffixes= c('.M', '.Tot'))

keep<- rowSums(cmp[,list(ear043_M8oxBS.Tot, ear045_T3oxBS.Tot)]) > 40
cmp<- cmp[keep, ]

## PRIOR: BETA BIMODAL
## ===================
# Estmate prior beta params for pct_met > x and pct_met < x by matching quantiles
cutoff<- 0.2
priorBetaParam_T_High<- fitdist(bdg[library_id == 'ear045_T3oxBS' & pct_met > cutoff & !is.na(pct_met), pct_met], 'beta', method= 'qme', probs= c(0.1, 0.9))
priorBetaParam_T_Low<- fitdist(bdg[library_id == 'ear045_T3oxBS' & pct_met < cutoff & !is.na(pct_met), pct_met], 'beta', method= 'qme', probs= c(0.1, 0.9))
priorBetaParam_M_High<- fitdist(bdg[library_id == 'ear043_M8oxBS' & pct_met > cutoff & !is.na(pct_met), pct_met], 'beta', method= 'qme', probs= c(0.1, 0.9))
priorBetaParam_M_Low<- fitdist(bdg[library_id == 'ear043_M8oxBS' & pct_met < cutoff & !is.na(pct_met), pct_met], 'beta', method= 'qme', probs= c(0.1, 0.9))

# Plot priors
prop_T_high<- smryH[library_id == 'ear045_T3oxBS', cnt] / (smryH[library_id == 'ear045_T3oxBS', cnt] + smryL[library_id == 'ear045_T3oxBS', cnt])
prop_M_high<- smryH[library_id == 'ear043_M8oxBS', cnt] / (smryH[library_id == 'ear043_M8oxBS', cnt] + smryL[library_id == 'ear043_M8oxBS', cnt])

rndPrior_T<- c(
    rbeta(10 * NSIM * prop_T_high, priorBetaParam_T_High$estimate[1], priorBetaParam_T_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_T_high), priorBetaParam_T_Low$estimate[1], priorBetaParam_T_Low$estimate[2]))
rndPrior_M<- c(
    rbeta(10 * NSIM * prop_M_high, priorBetaParam_M_High$estimate[1], priorBetaParam_M_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_M_high), priorBetaParam_M_Low$estimate[1], priorBetaParam_M_Low$estimate[2]))

pdf('bimodalBetaPrior.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(mgp= c(1.75, 0.5, 0), las= 1)
plot(density(rndPrior_T), xlim= c(0, 1), col= 'red', lwd= 2, xlab= '% methylated', main= 'Prior distribution')
points(density(rndPrior_M), col= 'blue', type= 'l', lwd= 2)
grid(col= 'grey30')
legend('topleft', lwd= 2, col= c('red', 'blue'), legend= c('Tumor', 'Margin'))
dev.off()
system('scp bimodalBetaPrior.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151030_betaBinEstimate/')

## DATA: Uniform prior
## ===================
yT = x$ear045_T3oxBS.M   # Tumor C
nT = x$ear045_T3oxBS.Tot # Tumor C+T

yM = x$ear043_M8oxBS.M   # Margin C
nM = x$ear043_M8oxBS.Tot # Margin C+T

data_T<- c(
    rbeta(NSIM * prop_T_high, yT + 1, (nT - yT) + 1),
    rbeta(NSIM * (1-prop_T_high), yT + 1, (nT - yT) + 1))
data_M<- c(
    rbeta(NSIM * prop_M_high, yM + 1, (nM - yM) + 1),
    rbeta(NSIM * (1-prop_M_high), yM + 1, (nM - yM) + 1))

pdf('example.tmp.pdf', w= 16/2.54, h= 12/2.54, pointsize= 10)
par(mgp= c(1.75, 0.5, 0), las= 0)
hist(data_M, breaks= 30, border= 'white', col= '#0000FF90', xlim= c(0, 1), freq= FALSE)
hist(data_T, breaks= 30, col= '#FF000090', border= 'white', add= TRUE, freq= FALSE)
abline(v= yT/nT, col= 'red', lty= 'dashed', lwd= 2)
abline(v= yM/nM, col= 'blue', lty= 'dashed', lwd= 2)
legend('topleft', lwd= 2, col= c('red', 'blue'), legend= c('Tumor', 'Margin'))
grid(col= 'grey30')
dev.off()
system(sprintf('scp example.tmp.pdf $mac_office:/media/groups/Research/sblab/public_folders/berald01/projects/20150501_methylation_brain/20151030_betaBinEstimate/exampleDataObs.%s-%s.%s-%s.pdf', yT, nT, yM, nM))

## POSTERIOR: Update high and low priors
## =====================================
post_T<- c(
    rbeta(NSIM * prop_T_high, yT + priorBetaParam_T_High$estimate[1], (nT - yT) + priorBetaParam_T_High$estimate[2]),
    rbeta(NSIM * (1-prop_T_high), yT + priorBetaParam_T_Low$estimate[1], (nT - yT) + priorBetaParam_T_Low$estimate[2]))
post_M<- c(
    rbeta(NSIM * prop_M_high, yM + priorBetaParam_M_High$estimate[1], (nM - yM) + priorBetaParam_M_High$estimate[2]),
    rbeta(NSIM * (1-prop_M_high), yM + priorBetaParam_M_Low$estimate[1], (nM - yM) + priorBetaParam_M_Low$estimate[2]))

pdf('example.tmp.pdf', w= 16/2.54, h= 12/2.54, pointsize= 10)
par(mgp= c(1.75, 0.5, 0), las= 0)
hist(post_M, breaks= 30, border= 'white', col= '#0000FF90', xlim= c(0, 1), freq= FALSE)
hist(post_T, breaks= 30, col= '#FF000090', border= 'white', add= TRUE, freq= FALSE)
abline(v= yT/nT, col= 'red', lty= 'dashed', lwd= 2)
abline(v= yM/nM, col= 'blue', lty= 'dashed', lwd= 2)
grid(col= 'grey30')
legend('topleft', lwd= 2, col= c('red', 'blue'), legend= c('Tumor', 'Margin'))
dev.off()
system(sprintf('scp example.tmp.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151030_betaBinEstimate/examplePostObs.%s-%s.%s-%s.pdf', yT, nT, yM, nM))

## DIFFERENCE
pdf('example.tmp.pdf', w= 16/2.54, h= 12/2.54, pointsize= 10)
par(mgp= c(1.75, 0.5, 0), las= 0)
hist(post_T - post_M, breaks= 30, border= 'white', col= 'grey60', xlim= c(-1, 1), freq= FALSE)
abline(v= 0, col= 'black', lty= 'dashed', lwd= 1)
grid(col= 'grey30')
legend('topleft', lwd= 2, col= c('red', 'blue'), legend= c('Tumor', 'Margin'))
dev.off()
system(sprintf('scp example.tmp.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151030_betaBinEstimate/exampleDiffObs.%s-%s.%s-%s.pdf', yT, nT, yM, nM))

## Posterior for all sites
## =======================
p<- c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995)
simPostDiff<- function(x, priorBetaParam_T_High, priorBetaParam_T_Low, priorBetaParam_M_High, priorBetaParam_M_Low,
    prop_T_high, prop_M_high){
    NSIM<- 100000
    x<- unlist(x)

    yM = x[1]   # Met in Margin
    yT = x[2]   # Met in tumor
    nM = x[3]   # Counts
    nT = x[4]

    post_T<- c(
        rbeta(NSIM * prop_T_high, yT + priorBetaParam_T_High$estimate[1], (nT - yT) + priorBetaParam_T_High$estimate[2]),
        rbeta(NSIM * (1-prop_T_high), yT + priorBetaParam_T_Low$estimate[1], (nT - yT) + priorBetaParam_T_Low$estimate[2]))
    post_M<- c(
        rbeta(NSIM * prop_M_high, yM + priorBetaParam_M_High$estimate[1], (nM - yM) + priorBetaParam_M_High$estimate[2]),
        rbeta(NSIM * (1-prop_M_high), yM + priorBetaParam_M_Low$estimate[1], (nM - yM) + priorBetaParam_M_Low$estimate[2]))

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

datOxBS<- cmp[, list(ear043_M8oxBS.M, ear045_T3oxBS.M, ear043_M8oxBS.Tot, ear045_T3oxBS.Tot)] ## Order of columns matters!
clus<- makeCluster(24)
clusterExport(clus, list('simPostDiff', 'Mode', 'priorBetaParam_T_High', 'priorBetaParam_T_Low', 'priorBetaParam_M_High', 'priorBetaParam_M_Low',
    'prop_T_high', 'prop_M_high', 'p'))
outqq<- parRapply(clus, datOxBS, function(x, y, z) simPostDiff(x, priorBetaParam_T_High, priorBetaParam_T_Low, priorBetaParam_M_High, priorBetaParam_M_Low,
    prop_T_high, prop_M_high))
stopCluster(clus)
save(outqq, file= 'outqq.tmp.Rdata')

# NB: It's not very smart to go thorugh the long list outqq 3 times!
postDT<- data.table(do.call('rbind', lapply(outqq, function(x) c(x[['post_T']], library_id= 'ear045_T3oxBS')))) # 5mC in tumor
postDT<- rbindlist(list(data.table(do.call('rbind', lapply(outqq, function(x) c(x[['post_M']], library_id= 'ear043_M8oxBS')))), postDT)) # 5mC in margin
postDT<- rbindlist(list(data.table(do.call('rbind', lapply(outqq, function(x) c(x[['post_D']], library_id= 'oxBS_TvsM')))), postDT)) # 5mC difference T-M
for(x in p){
    pp<- paste('p', x, sep= '')
    postDT[[pp]]<- as.numeric(postDT[[pp]])
}
postDT[['mode']]<- as.numeric(postDT[['mode']])
stopifnot(nrow(postDT) == nrow(cmp) * 3)
postDT[, chrom := rep(cmp$chrom, 3)]
postDT[, start := rep(cmp$start, 3)]
postDT[, end := rep(cmp$end, 3)]

oxBS_post<- postDT[library_id %in% c('ear043_M8oxBS', 'ear045_T3oxBS'),]
write.table(file= 'oxBS_posterior.txt', x= oxBS_post, quote= FALSE, row.names= FALSE, sep= '\t')
system('pigz -f oxBS_posterior.txt')

oxBS_TvsM<- postDT[library_id == 'oxBS_TvsM',]
oxBS_TvsM[, library_id := NULL]

## Infer the quantile in the posterior distr where the 0 point is crossed
## This part might be included in the parallel loop!!
## ======================================================================
getQAtZero<- function(x, p){
    xlm<- glm(p ~ x, family= quasibinomial)
    qz<- inv.logit(xlm$coefficients[1])
    return(qz)
}
dat<- as.matrix((oxBS_TvsM[, paste('p', p, sep= ''), with= FALSE]))
clus<- makeCluster(24)
clusterExport(clus, list('getQAtZero', 'inv.logit', 'p'))
q0<- parRapply(clus, dat, function(x) getQAtZero(x, p))
stopCluster(clus)

## Convert quantile to a p-value stats. The smallest the more far in the tail]
## the zero point is.
oxBS_TvsM[, pqAt0 := ifelse(q0 > 0.5, 1-q0, q0)]

write.table(file= 'oxBS_TvsM_posterior.txt', x= oxBS_TvsM, quote= FALSE, row.names= FALSE, sep= '\t')
system('pigz -f oxBS_TvsM_posterior.txt')
```

## Summary plots

```R
## Hist T - M difference
## =====================
gg<- ggplot(data= oxBS_TvsM, aes(x= mode)) +
    geom_histogram(colour= 'white', fill= 'blue', alpha= 0.3) +
    geom_histogram(data= oxBS_TvsM[pqAt0 < 0.01], aes(x= mode), colour= 'white', fill= 'red', alpha= 0.3) +
    xlab('5mC [Tumor - Margin]\nPosterior mode') +
    ggtitle('Difference in 5mC between tumor and margin') +
    scale_x_continuous(breaks=pretty_breaks(n=10)) +
    geom_vline(xintercept= 0, colour= 'grey30', linewidth= 2, linetype= 'dotted')
ggsave('5mC_TvsM_hist.pdf', gg, width= 16, height= 12, units= 'cm')
system('rsync --remove-source-files 5mC_TvsM_hist.pdf $mac_office:$cri_public_projects/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/20151102_Tumor_vs_Margin_by_bayesian_approach/')

## Plot width of CI vs mean depth
## ==============================
oxBS_TvsM<- merge(oxBS_TvsM, cmp, by= c('chrom', 'start', 'end'))
xat<- seq(1, nrow(oxBS_TvsM), length.out= 10000)
gg<- ggplot(data= oxBS_TvsM[xat,], aes(x= rowMeans(oxBS_TvsM[xat, list(ear043_M8oxBS.Tot, ear045_T3oxBS.Tot)]), y= oxBS_TvsM[xat, p0.995 - p0.005])) +
    geom_point(alpha= 0.2, size= 1) +
    geom_smooth(se= FALSE) +
    xlim(NA, 100) +
    xlab('Depth [mean(tumor, margin)]') +
    ylab('Width CI') +
    ggtitle('Difference in Tumor - Margin 5mC\nWidth of credible interval in relation to sequencing depth')
ggsave('5mC_TvsM_width_CI.pdf', width= 12, height= 12, units= 'cm')
system('rsync --remove-source-files 5mC_TvsM_width_CI.pdf $mac_office:$cri_public_projects/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/20151102_Tumor_vs_Margin_by_bayesian_approach/')

## Plot posterior T and M vs observed
## ==================================
oxBS_post<- merge(oxBS_post, bdg, by= c('chrom', 'start', 'end', 'library_id'))
xat<- seq(1, nrow(oxBS_post), length.out= 50000)
gg<- ggplot(data= oxBS_post[xat], aes(x= 100 * pct_met, y= 100 * (pct_met - p0.5), colour= cnt_tot)) +
    geom_point(alpha= 0.2, size= 1) +
    facet_wrap(~library_id) +
    geom_hline(yintercept= 0, colour= 'blue', alpha= 0.2, linewidth= 2) +
    scale_colour_gradient2(limits=c(5, 60), midpoint= 35) + # low="blue", mid= 'green', high= 'red', 
    xlab('Observed % 5mC') +
    ylab('% 5mC [Obs - Posterior]') +
    ggtitle('Posterior estimate of 5mC vs observed')
ggsave('5mC_posterior_vs_obs.pdf', width= 20, height= 10, units= 'cm')
system('rsync --remove-source-files 5mC_posterior_vs_obs.pdf $mac_office:$cri_public_projects/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/20151102_Tumor_vs_Margin_by_bayesian_approach/')

## Tumor vs Margin 5mC
## ===================
dat<- dcast.data.table(data= bdg, start ~ library_id, value.var= 'pct_met')
xat<- seq(1, nrow(dat), length.out= 10000)
ggRaw<- ggplot(data= dat[xat,], aes(100 * ear043_M8oxBS, y= 100 * (ear045_T3oxBS - ear043_M8oxBS))) +
    geom_point(alpha= 0.2, size= 1) +
    geom_hline(yintercept= 0, colour= 'blue', alpha= 0.2, linewidth= 2) +
    xlim(0, 100) +
    ylim(-100, 100) +
    xlab('Margin %5mC') +
    ylab('% 5mC [Tumor - Margin]') +
    ggtitle('Margin vs tumor 5mC - Observed')
dat<- dcast.data.table(data= oxBS_post, start ~ library_id, value.var= 'mode')
xat<- seq(1, nrow(dat), length.out= 10000)
ggPost<- ggplot(data= dat[xat,], aes(100 * ear043_M8oxBS, y= 100 * (ear045_T3oxBS - ear043_M8oxBS))) +
    geom_point(alpha= 0.2, size= 1) +
    geom_hline(yintercept= 0, colour= 'blue', alpha= 0.2, linewidth= 2) +
    xlim(0, 100) + 
    ylim(-100, 100) +
    xlab('Margin %5mC') +
    ylab('% 5mC [Tumor - Margin]') +
    ggtitle('Margin vs tumor 5mC - Posterior estimates')
gg<- arrangeGrob(ggRaw, ggPost, nrow= 2)
ggsave('5mC_TvsM.pdf', gg, width= 16, height= 18, units= 'cm')
system('rsync --remove-source-files 5mC_TvsM.pdf $mac_office:$cri_public_projects/20150921_BrainMethylomeRoadMap/20151013_methylation_chr18/20151102_Tumor_vs_Margin_by_bayesian_approach/')
```

