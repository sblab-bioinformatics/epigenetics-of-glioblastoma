<!-- MarkdownTOC -->

- [Estimating methylation and hydroxymethylation at base resolution](#estimating-methylation-and-hydroxymethylation-at-base-resolution)
- [Bayesian estimate of 5mC in tumor and margin](#bayesian-estimate-of-5mc-in-tumor-and-margin)
- [Bayesian estimate of differential methylation](#bayesian-estimate-of-differential-methylation)

<!-- /MarkdownTOC -->

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