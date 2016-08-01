CUTOFF<- 0.4
pct_met<- bdg[cnt_met_mar/cnt_tot_mar > CUTOFF, cnt_met_mar/cnt_tot_mar]
priorBetaParam_M_High<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

pct_met<- bdg[cnt_met_mar/cnt_tot_mar <= CUTOFF, cnt_met_mar/cnt_tot_mar]
priorBetaParam_M_Low<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

## Tumor, high and low
pct_met<- bdg[cnt_met_tum/cnt_tot_tum > CUTOFF, cnt_met_tum/cnt_tot_tum]
priorBetaParam_T_High<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 
    'beta', method= 'qme', probs= c(0.1, 0.9))

pct_met<- bdg[cnt_met_tum/cnt_tot_tum <= CUTOFF, cnt_met_tum/cnt_tot_tum]
priorBetaParam_T_Low<- fitdist(pct_met[sample(1:length(pct_met), 
    size= length(pct_met)/100)], 'beta', method= 'qme', probs= c(0.1, 0.9))

prop_M_high<- nrow(bdg[cnt_met_mar/cnt_tot_mar > CUTOFF]) / nrow(bdg)
prop_T_high<- nrow(bdg[cnt_met_tum/cnt_tot_tum > CUTOFF]) / nrow(bdg)

NSIM<- 1000000
rndPrior_M<- c(
    rbeta(10 * NSIM * prop_M_high, priorBetaParam_M_High$estimate[1], priorBetaParam_M_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_M_high), priorBetaParam_M_Low$estimate[1], priorBetaParam_M_Low$estimate[2]))
rndPrior_T<- c(
    rbeta(10 * NSIM * prop_T_high, priorBetaParam_T_High$estimate[1], priorBetaParam_T_High$estimate[2]),
    rbeta(10 * NSIM * (1-prop_T_high), priorBetaParam_T_Low$estimate[1], priorBetaParam_T_Low$estimate[2]))

betaPars<- t(data.frame(
        margin_high= priorBetaParam_M_High$estimate,
        margin_low= priorBetaParam_M_Low$estimate,
        tumor_high= priorBetaParam_T_High$estimate,
        tumor_low= priorBetaParam_T_Low$estimate))
print(betaPars)

Mode<- function(x) {
  ## Helper function to get mode of vector x. I.e. find the max of the 
  ## kernel density of x
  ## Default bandwidth is not the best, but it's fast
  z<- density(x)    
  return( z$x[z$y==max(z$y)] )
}


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

    ## Vector of quantile probabilities.
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

    ## Quantile and mode of tumor posterior
    Tq<- quantile(post_T, p);
    names(Tq)<- nn
    Tq<- c(Tq, mode= Mode(post_T), avg= mean(post_T))

    ## Quantile and mode of margin posterior
    Mq<- quantile(post_M, p)
    names(Mq)<- nn
    Mq<- c(Mq, mode= Mode(post_M), avg= mean(post_M))

    ## Quantiles and mode of the difference
    Dq<- quantile(post_T - post_M, p)
    names(Dq)<- nn
    Dq<- c(Dq, mode= Mode(post_T - post_M), avg= mean(post_T - post_M))

    qq<- list(post_T= Tq, post_M= Mq, post_D= Dq)
    return(qq)
}

## Apply simPostDiff to all CpGs:
xsub<- seq(1, nrow(bdg), length.out= 2000)
datOxBS<- bdg[xsub, list(cnt_met_mar,
                     cnt_tot_mar,
                     cnt_met_tum,
                     cnt_tot_tum)] ## Order of columns matters!
clus<- makeCluster(8)
clusterExport(clus, list('simPostDiff', 'Mode', 'priorBetaParam_T_High', 'priorBetaParam_T_Low', 'priorBetaParam_M_High', 'priorBetaParam_M_Low',
    'prop_T_high', 'prop_M_high'))
system.time({
outqq<- parRapply(clus, datOxBS, 
    function(x) simPostDiff(x, 
        priorBetaParam_T_High, priorBetaParam_T_Low, priorBetaParam_M_High, priorBetaParam_M_Low,
        prop_T_high, prop_M_high, nsim= 200000))
})
stopCluster(clus)

## Convert to data.table
cnames<- names(outqq[[1]][[1]])
slots<- names(outqq[[1]])
postDT<- data.table(matrix(unlist(outqq), ncol= length(cnames) * length(slots), byrow= TRUE))
rm(outqq); gc()

xnames<-  apply(expand.grid(cnames, slots), 1, paste, collapse= '_')
setnames(postDT, names(postDT), xnames)

# Add to posterior datatable the position and raw counts
postDT<- cbind(bdg[xsub], postDT)
gc()