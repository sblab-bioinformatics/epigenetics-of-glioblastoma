library(scales)
library(ggplot2)
library(data.table)
library(parallel)
library(fitdistrplus)

getFisherPval<- function(x){
    # x: Named list to get counts for fisher test. Note that you input total 
    # counts which will be converted to count unmethylated for fisher test.
    # E.g.
    # x<- list(cnt_met_mar= 14, cnt_tot_mar= 22, cnt_met_tum= 8, cnt_tot_tum= 11)
    mat<- matrix(data= 
        c(x[['cnt_met_mar']], x[['cnt_tot_mar']] - x[['cnt_met_mar']],
          x[['cnt_met_tum']], x[['cnt_tot_tum']] - x[['cnt_met_tum']]),
        nrow= 2)
    return(fisher.test(mat, conf.int= FALSE)$p.value)
}

getBetaParams<- function(mu, sigma){
    # Return the alpha and beta params for a beta distribution with desired mean 
    # and stdev. Given mean and stdev of a beta distr, the corresponding alpha 
    # and beta params are:
    # ⎧      ⎛ 2        2⎞              ⎛ 2        2⎞⎫
    # ⎪   -μ⋅⎝μ  - μ + σ ⎠      (μ - 1)⋅⎝μ  - μ + σ ⎠⎪
    # ⎨α: ─────────────────, β: ─────────────────────⎬
    # ⎪            2                       2         ⎪
    # ⎩           σ                       σ          ⎭
    stopifnot(sigma > 0)    
    alpha<- (-mu * (mu^2 - mu + sigma^2)) / sigma^2
    beta<- ( (mu - 1) * (mu^2 - mu + sigma^2) ) / sigma^2
    stopifnot(alpha > 0)
    stopifnot(beta > 0)
    return(list(a= alpha, b= beta))
}

inv.logit<- function(x){
    exp(x)/(1+exp(x))
}
getQAtZero<- function(x, p){
    # Fit a (quasi)binomial regression as prob ~ quantile. Then extract the
    # probability where the quantile is zero.
    xlm<- glm(p ~ x, family= quasibinomial)
    qz<- inv.logit(xlm$coefficients[1])
    return(qz)
}
