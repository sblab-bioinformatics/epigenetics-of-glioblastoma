<!-- MarkdownTOC -->

- [A Bayesian estimate of the 5hmC bsgin](#a-bayesian-estimate-of-the-5hmc-bsgin)
    - [Prepare data](#prepare-data)
    - [Prepare prior distribution](#prepare-prior-distribution)
    - [Posterior 5mC \(oxBS\) and 5hmC \(BS - oxBS\)](#posterior-5mc-oxbs-and-5hmc-bs---oxbs)
    - [Convert quantile to a p-value stats. The smallest the more far in the tail](#convert-quantile-to-a-p-value-stats-the-smallest-the-more-far-in-the-tail)
    - [the zero point is.](#the-zero-point-is)
    - [NB: Previous file overwritten!](#nb-previous-file-overwritten)

<!-- /MarkdownTOC -->


A Bayesian estimate of the 5hmC bsgin
======================================

Closely following the method in [bayesian_estimate_of_5mC_5hmC](bayesian_estimate_of_5mC_5hmC.md)


Prepare data
------------

<!---
cd /nas/sblab_data1/group_folders/berald01/projects/20150501_methylation_brain/20160613_bayes/
-->


```R
R
library(ggplot2)
library(data.table)
library(gridExtra)
library(parallel)
library(fitdistrplus)

## Data preparation
bs<- fread('zcat ear042_M8BS.cpg.bedGraph.gz') 
oxbs<- fread('zcat ear043_M8oxBS.cpg.bedGraph.gz')

stopifnot(all(bs$V3 - bs$V2 == 2))
stopifnot(all(oxbs$V3 - oxbs$V2 == 2))

bs[, V3 := NULL]
oxbs[, V3 := NULL]

xn<- c('chrom', 'start', 'cnt_met', 'cnt_tot')
setnames(bs, names(bs), xn)
setnames(oxbs, names(oxbs), xn)
bdg<- merge(bs, oxbs, by= c('chrom', 'start'), suffixes= c('_bs', '_oxbs'))
rm(bs, oxbs)
```

Methylation tends to follow a bimodal distribution:

<!--
oxBS -> mar \
             } Difference expected +ve
BS<- tum    /
-->

```R
set.seed(1234)
xn<- sample(1:nrow(bdg), size= nrow(bdg)/100, replace= FALSE)
gg<- ggplot(data= bdg[xn]) +
    geom_histogram(aes(x= 100 * cnt_met_oxbs/cnt_tot_oxbs), fill= 'blue', alpha= 0.4) +
    geom_histogram(aes(x= 100 * cnt_met_bs/cnt_tot_bs), fill= 'red', alpha= 0.4) +
    ylab('N. CpG x1000') +
    xlab('% 5mC') +
    annotate("text", x = c(0, 0), y = c(40000, 38000), label = c("oxBS", 'BS'), 
        colour= c('blue', 'red'), vjust= 1) +
    ggtitle('Distribution of C-modification in margin')
ggsave('hist_bs_oxbs_margin.png', w= 14, h= 12, units= 'cm')
```

Prepare prior distribution
--------------------------


```R
CUTOFF<- 0.5

## oxBS, high and low
pct_met<- bdg[cnt_met_oxbs/cnt_tot_oxbs > CUTOFF, cnt_met_oxbs/cnt_tot_oxbs]
priorBetaParam_oxbs_High<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

pct_met<- bdg[cnt_met_oxbs/cnt_tot_oxbs <= CUTOFF, cnt_met_oxbs/cnt_tot_oxbs]
priorBetaParam_oxbs_Low<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

## BS, high and low
pct_met<- bdg[cnt_met_bs/cnt_tot_bs > CUTOFF, cnt_met_bs/cnt_tot_bs]
priorBetaParam_bs_High<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

pct_met<- bdg[cnt_met_bs/cnt_tot_bs <= CUTOFF, cnt_met_bs/cnt_tot_bs]
priorBetaParam_bs_Low<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))
```

       Parameter | oxBS   | BS 
---------------- | ------ | -----
α<sub>low</sub>  | 0.25   | 0.21
β<sub>low</sub>  | 1.48   | 1.42
α<sub>high</sub> | 11.72  | 11.74
β<sub>high</sub> | 4.97   | 2.55


Plot prior

```R
prop_oxbs_high<- nrow(bdg[cnt_met_oxbs/cnt_tot_oxbs > CUTOFF]) / nrow(bdg)
prop_bs_high<- nrow(bdg[cnt_met_bs/cnt_tot_bs > CUTOFF]) / nrow(bdg)

NSIM<- 100000
rndPrior_oxbs<- c(
    rbeta(10 * NSIM * prop_oxbs_high, priorBetaParam_oxbs_High$estimate[1], priorBetaParam_oxbs_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_oxbs_high), priorBetaParam_oxbs_Low$estimate[1], priorBetaParam_oxbs_Low$estimate[2]))
rndPrior_bs<- c(
    rbeta(10 * NSIM * prop_bs_high, priorBetaParam_bs_High$estimate[1], priorBetaParam_bs_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_bs_high), priorBetaParam_bs_Low$estimate[1], priorBetaParam_bs_Low$estimate[2]))

## Plot priors:
gg<- ggplot(data= NULL) +
    geom_histogram(aes(x= 100 * rndPrior_oxbs, y= ..density..), fill= 'blue', alpha= 0.4) +
    geom_histogram(aes(x= 100 * rndPrior_bs, y= ..density..), fill= 'red', alpha= 0.4) +
    xlab('% 5mC') + 
    ggtitle('Prior distributions of 5mC in oxbsgin and bsor')
ggsave('bimodalBetaPrior_bs_oxbs_margin.png', w= 14, h= 12, units= 'cm')
```


Posterior 5mC (oxBS) and 5hmC (BS - oxBS)
----------------------------------------

The posterior 5mC is redundant as it is the same as in [bayesian_estimate_of_5mC_5hmC](bayesian_estimate_of_5mC_5hmC.md)

```R
Mode<- function(x) {
  ## Helper function to get mode of vector x. I.e. find the max of the 
  ## kernel density of x
  ## Default bandwidth is not the best, but it's fast
  z<- density(x)    
  return( z$x[z$y==max(z$y)] )
}


simPostDiff<- function(x, priorBetaParam_bs_High, priorBetaParam_bs_Low, 
                          priorBetaParam_oxbs_High, priorBetaParam_oxbs_Low,
                          prop_bs_high, prop_oxbs_high, nsim= 20000){
    ## Update priors given data in x
    ## x: 
    ##      Vectors of length 4 giving the observed data: 
    ##      1) count methylated in oxbs 2) count total in oxbs
    ##      3) count methylated bs 4) count total bs
    ## priorBetaParam_bs/Low:
    ##      Parameters of the beta distribution of the low and high mC tails, 
    ##      for bs and oxbs
    ## prop_bs/high: 
    ##      Proportion of the entire mC distribution belonging to the 'high' tail
    ## nsim:
    ##      Number of random samples to draw for simulation
    ## Returns:
    ##      List of 3 vectors: Quantiles and mode of the posteriors of: 
    ##      1) bs 5mC, 2) oxbs 5mC, 3) difference BS-oxBS

    ## Vector of quantile probabilities.
    p<- c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995)
    
    x<- unlist(x)
    stopifnot(length(x) == 4)
    stopifnot(all(!is.na(x)))

    yM = x[1]   # Met in oxbs
    nM = x[2]   # Counts
    yT = x[3]   # Met in bs
    nT = x[4]

    ## Some simple sanity check
    stopifnot(yT <= nT)
    stopifnot(yM <= nM)

    post_bs<- c(
        rbeta(nsim * prop_bs_high, yT + priorBetaParam_bs_High$estimate[1], (nT - yT) + priorBetaParam_bs_High$estimate[2]),
        rbeta(nsim * (1-prop_bs_high), yT + priorBetaParam_bs_Low$estimate[1], (nT - yT) + priorBetaParam_bs_Low$estimate[2]))
    post_oxbs<- c(
        rbeta(nsim * prop_oxbs_high, yM + priorBetaParam_oxbs_High$estimate[1], (nM - yM) + priorBetaParam_oxbs_High$estimate[2]),
        rbeta(nsim * (1-prop_oxbs_high), yM + priorBetaParam_oxbs_Low$estimate[1], (nM - yM) + priorBetaParam_oxbs_Low$estimate[2]))

    nn<- sprintf('p%s', p)
    
    ## Quantile and mode of bsor posterior
    Tq<- quantile(post_bs, p);
    names(Tq)<- nn
    Tq<- c(Tq, mode= Mode(post_bs))
    
    ## Quantile and mode of oxbsgin posterior
    Mq<- quantile(post_oxbs, p)
    names(Mq)<- nn
    Mq<- c(Mq, mode= Mode(post_oxbs))
    
    ## Quantiles and mode of the difference
    Dq<- quantile(post_bs - post_oxbs, p)
    names(Dq)<- nn
    Dq<- c(Dq, mode= Mode(post_bs - post_oxbs))
    
    qq<- list(post_bs= Tq, post_oxbs= Mq, post_D= Dq)
    return(qq)
}

## Apply simPostDiff to all CpGs:
datOxBS<- bdg[, list(cnt_met_oxbs,
                     cnt_tot_oxbs,
                     cnt_met_bs,
                     cnt_tot_bs)] ## Order of columns matters!
clus<- makeCluster(24)
clusterExport(clus, list('simPostDiff', 'Mode', 'priorBetaParam_bs_High', 'priorBetaParam_bs_Low', 'priorBetaParam_oxbs_High', 'priorBetaParam_oxbs_Low',
    'prop_bs_high', 'prop_oxbs_high'))
system.time({
outqq<- parRapply(clus, datOxBS, 
    function(x) simPostDiff(x, 
        priorBetaParam_bs_High, priorBetaParam_bs_Low, priorBetaParam_oxbs_High, priorBetaParam_oxbs_Low,
        prop_bs_high, prop_oxbs_high, nsim= 20000))
})
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

# You might want to quit R and reload this obj to clean up the memory. You need to reload also bdg
postDT<- fread('postDT.tmp.txt')

# Add to posterior datatable the position and raw counts
postDT<- cbind(bdg, postDT)
rm(bdg); gc()

# Write out tables
write.table(x= 
  postDT[, list(chrom,
                start,
                cnt_met_oxbs,
                cnt_tot_oxbs,
                cnt_met_bs,
                cnt_tot_bs,
                p0.005= 100 * round(p0.005_post_D, 4),
                p0.025= 100 * round(p0.025_post_D, 4),
                p0.25=  100 * round(p0.25_post_D, 4),
                p0.5=   100 * round(p0.5_post_D, 4),
                p0.75=  100 * round(p0.75_post_D, 4),
                p0.975= 100 * round(p0.975_post_D, 4),
                p0.995= 100 * round(p0.995_post_D, 4),
                mode= 100 * round(mode_post_D, 4))],
  file= 'posterior_5hmC_margin.txt', row.names= FALSE, quote= FALSE, sep= '\t')

write.table(x= 
  postDT[, list(chrom,
                start,
                p0.005= 100 * round(p0.005_post_bs, 4),
                p0.025= 100 * round(p0.025_post_bs, 4),
                p0.25=  100 * round(p0.25_post_bs, 4),
                p0.5=   100 * round(p0.5_post_bs, 4),
                p0.75=  100 * round(p0.75_post_bs, 4),
                p0.975= 100 * round(p0.975_post_bs, 4),
                p0.995= 100 * round(p0.995_post_bs, 4),
                mode= 100 * round(mode_post_bs, 4))],
  file= 'posterior_bs_margin.txt', row.names= FALSE, quote= FALSE, sep= '\t')

write.table(x= 
  postDT[, list(chrom,
                start,
                p0.005= 100 * round(p0.005_post_oxbs, 4),
                p0.025= 100 * round(p0.025_post_oxbs, 4),
                p0.25=  100 * round(p0.25_post_oxbs, 4),
                p0.5=   100 * round(p0.5_post_oxbs, 4),
                p0.75=  100 * round(p0.75_post_oxbs, 4),
                p0.975= 100 * round(p0.975_post_oxbs, 4),
                p0.995= 100 * round(p0.995_post_oxbs, 4),
                mode= 100 * round(mode_post_oxbs, 4))],
  file= 'posterior_oxbs_margin.txt', row.names= FALSE, quote= FALSE, sep= '\t')
system('gzip posterior_*.txt')
```

Probability of the posterior distribution to include zero (*i.e.* probability that the difference between is not different from zero). 

```
<!--
This part should be included in the calculation of the posterior above!
-->
```R
R
library(data.table)
library(scales)
library(ggplot2)
library(parallel)

xfile<- 'posterior_5hmC_margin.txt'
postDiff<- fread(sprintf('zcat %s.gz', xfile))

inv.logit<- function(x){
    exp(x)/(1+exp(x))
}
getQAtZero<- function(x, p){
    ## Fit a (quasi)binomial regression as prob ~ quantile. 
    ## Then extract the probability where the quantile is zero.
    xlm<- glm(p ~ x, family= quasibinomial)
    qz<- inv.logit(xlm$coefficients[1])
    return(qz)
}
dat<- as.matrix((postDiff[, list(p0.005, p0.025, p0.25, p0.5, p0.75, p0.975, p0.995)]))
p<- c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995)
clus<- makeCluster(24)
clusterExport(clus, list('getQAtZero', 'inv.logit', 'p'))
q0<- parRapply(clus, dat, function(x) getQAtZero(x, p))
stopCluster(clus)

## Convert quantile to a p-value stats. The smallest the more far in the tail
## the zero point is.
postDiff[, prob0 := ifelse(q0 > 0.5, 1-q0, q0)]

## NB: Previous file overwritten!
write.table(x= postDiff, file= xfile, row.names= FALSE, quote= FALSE, sep= '\t')
system(sprintf('gzip %s', xfile))
```

Inspection of 5hmC
------------------


```R
R
library(data.table)
library(scales)
library(ggplot2)

postDiff<- fread('zcat posterior_5hmC_margin.txt.gz')
xn<- seq(1, nrow(postDiff), length.out= 50000)
gg<- ggplot(data= postDiff[xn, ], aes(x= 100 * ((cnt_met_bs/cnt_tot_bs) - (cnt_met_oxbs/cnt_tot_oxbs)), y= mode,
    colour= ifelse(cnt_tot_oxbs + cnt_tot_bs > 200, 200, cnt_tot_oxbs + cnt_tot_bs))) +
    scale_colour_gradient2('Read count', low=muted("blue"), high=muted("red"), midpoint= 75, mid= muted('blue')) +
    geom_point(alpha= 0.50, size= 0.15) +
    geom_abline(intercept= 0, slope= 1, linetype= 'dashed') +
    xlab('Observed 5hmC difference') +
    ylab('Bayesian posterior 5hmC') +
    ggtitle('5hmC in margin\nObserved from raw counts vs Bayesian estimate')
ggsave('posterior_difference_h5mc_margin.png', w= 14, h= 12, units= 'cm')
```



