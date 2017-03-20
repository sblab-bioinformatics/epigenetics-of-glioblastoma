
library(data.table)
library(ggplot2)

nucleotide.freq.cpg.C01 <- fread("cat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_C01.txt")
setnames(nucleotide.freq.cpg.C01, c("chr", "start", "ref", "nA", "nC", "nG", "nT"))

nucleotide.freq.cpg.E01 <- fread("cat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_E01.txt")
setnames(nucleotide.freq.cpg.E01, c("chr", "start", "ref", "nA", "nC", "nG", "nT"))

ear042_M8BS.cpg <- fread("zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear042_M8BS.cpg.bedGraph.gz")
setnames(ear042_M8BS.cpg, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))
ear042_M8BS.cpg[, pct_met := cnt_met/cnt_tot]

ear043_M8oxBS.cpg <- fread("zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear043_M8oxBS.cpg.bedGraph.gz")
setnames(ear043_M8oxBS.cpg, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))
ear043_M8oxBS.cpg[, pct_met := cnt_met/cnt_tot]

ear044_T3BS.cpg <- fread("zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear044_T3BS.cpg.bedGraph.gz")
setnames(ear044_T3BS.cpg, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))
ear044_T3BS.cpg[, pct_met := cnt_met/cnt_tot]

ear045_T3oxBS.cpg <- fread("zcat /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear045_T3oxBS.cpg.bedGraph.gz")
setnames(ear045_T3oxBS.cpg, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand"))
ear045_T3oxBS.cpg[, pct_met := cnt_met/cnt_tot]

setkey(ear042_M8BS.cpg, chr, start, end, strand)
setkey(ear043_M8oxBS.cpg, chr, start, end, strand)
setkey(ear044_T3BS.cpg, chr, start, end, strand)
setkey(ear045_T3oxBS.cpg, chr, start, end, strand)

M8.cpg <- merge(ear042_M8BS.cpg, ear043_M8oxBS.cpg, suffixes = c(".BS", ".oxBS"))
M8.cpg[, pct_5hmC := pct_met.BS - pct_met.oxBS]
M8.cpg[, pct_5mC := pct_met.oxBS]
M8.cpg[, pct_C_other := 1 - pct_met.BS]

T3.cpg <- merge(ear044_T3BS.cpg, ear045_T3oxBS.cpg, suffixes = c(".BS", ".oxBS"))
T3.cpg[, pct_5hmC := pct_met.BS - pct_met.oxBS]
T3.cpg[, pct_5mC := pct_met.oxBS]
T3.cpg[, pct_C_other := 1 - pct_met.BS]

rm(ear042_M8BS.cpg, ear043_M8oxBS.cpg, ear044_T3BS.cpg, ear045_T3oxBS.cpg)

nucleotide.freq.cpg.C01[, "count_alternative" := ifelse(ref == "C", nA+nG+nT, nA+nC+nT)]
nucleotide.freq.cpg.C01[, "count_total" := nA+nC+nG+nT]
nucleotide.freq.cpg.C01[, "pct_alternative" := ifelse(ref == "C", (nA+nG+nT)/(nA+nC+nG+nT), (nA+nC+nT)/(nA+nC+nG+nT)), with=FALSE]
nucleotide.freq.cpg.E01[, "count_alternative" := ifelse(ref == "C", nA+nG+nT, nA+nC+nT)]
nucleotide.freq.cpg.E01[, "count_total" := nA+nC+nG+nT]
nucleotide.freq.cpg.E01[, "pct_alternative" := ifelse(ref == "C", (nA+nG+nT)/(nA+nC+nG+nT), (nA+nC+nT)/(nA+nC+nG+nT)), with=FALSE]

setkey(nucleotide.freq.cpg.C01, chr, start, ref)
setkey(nucleotide.freq.cpg.E01, chr, start, ref)
nucleotide.freq.cpg.C01.E01 <- merge(nucleotide.freq.cpg.C01, nucleotide.freq.cpg.E01, suffixes = c(".C01", ".E01"))
nucleotide.freq.cpg.C01.E01[, c("nA.C01", "nC.C01", "nG.C01", "nT.C01", "nA.E01", "nC.E01", "nG.E01", "nT.E01"):= NULL]
rm(nucleotide.freq.cpg.C01, nucleotide.freq.cpg.E01)

setkey(M8.cpg, chr, start, end, strand)
setkey(T3.cpg, chr, start, end, strand)
methylation.cpg.M8.T3 <- merge(M8.cpg, T3.cpg, suffixes = c(".M8", ".T3"))
rm(M8.cpg, T3.cpg)

nucleotide.freq.cpg.C01.E01[,chr:=paste("chr", chr, sep="")]
setkey(nucleotide.freq.cpg.C01.E01, chr, start)
setkey(methylation.cpg.M8.T3, chr, start)
nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3 <- merge(nucleotide.freq.cpg.C01.E01, methylation.cpg.M8.T3)

dim(nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3)

rm(nucleotide.freq.cpg.C01.E01, methylation.cpg.M8.T3)

CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph")
setnames(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic, c("chr", "start", "end", "ref", "alt", "qss", "tqss"))

setkey(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic, chr, start, end, ref)
setkey(nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3, chr, start, end, ref)
snv.nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3 <- merge(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic, nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3, all=TRUE)
dim(snv.nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3)

rm(CancerLP2000729_DNA_E01_NormalLP2000729_DNA_C01.somatic, nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3)

l = list(snv.nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3[alt!="NA" & pct_alternative.C01!="NA"], snv.nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3[is.na(alt),])
setattr(l, 'names', c("SNV", "non-SNV"))
snv.nonsnv <- rbindlist(l, use.names=TRUE, idcol="ID")
rm(l)
snv.nonsnv[, c("ref", "alt", "qss", "tqss", "count_alternative.C01", "count_total.C01", "pct_alternative.C01", "count_alternative.E01", "count_total.E01", "pct_alternative.E01", "strand", "cnt_met.BS.M8", "cnt_tot.BS.M8", "pct_met.BS.M8", "cnt_met.oxBS.M8", "cnt_tot.oxBS.M8", "pct_met.oxBS.M8", "cnt_met.BS.T3", "cnt_tot.BS.T3", "pct_met.BS.T3", "cnt_met.oxBS.T3", "cnt_tot.oxBS.T3", "pct_met.oxBS.T3", "pct_5hmC.T3", "pct_5mC.T3", "pct_C_other.T3"):= NULL]
rm(snv.nucleotide.freq.cpg.C01.E01_methylation.cpg.M8.T3)
gc()
snv.nonsnv.m <- reshape2::melt(snv.nonsnv, id.vars = c("ID", "chr", "start", "end"), variable.name = "type_met", value.name = "pct_met")
rm(snv.nonsnv)

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/figures/20160514_combined.snv.nonsnv.cpg.M8.boxplot.pdf", width= 14/2.54, height= 10/2.54)
ggplot(snv.nonsnv.m, aes(x = factor(type_met, levels(factor(type_met))[c(2,1,3)]), y = 100*pct_met, fill = factor(ID, levels(factor(ID))[c(2,1)]))) + geom_boxplot(outlier.shape = NA, aes(colour = factor(type_met, levels(factor(type_met))[c(2,1,3)]))) + theme_bw() + theme(legend.title = element_blank()) + scale_fill_manual(values=c("#999999", "#E69F00")) + scale_colour_manual(values = c("turquoise3", "green", "red"), labels=c("5mC","5hmC","C")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/figures/20160514_combined.snv.nonsnv.cpg.M8.boxplot.nocolour.pdf", width= 14/2.54, height= 10/2.54)
ggplot(snv.nonsnv.m, aes(x = factor(type_met, levels(factor(type_met))[c(2,1,3)]), y = 100*pct_met, fill = factor(ID, levels(factor(ID))[c(2,1)]))) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(legend.title = element_blank()) + scale_fill_manual(values=c("#999999", "#E69F00")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/figures/20160514_combined.snv.nonsnv.cpg.M8.violin.pdf", width= 14/2.54, height= 10/2.54)
ggplot(snv.nonsnv.m, aes(x = factor(type_met, levels(factor(type_met))[c(2,1,3)]), y = 100*pct_met, fill = factor(ID, levels(factor(ID))[c(2,1)]))) + geom_violin(scale = "area", aes(colour = factor(type_met, levels(factor(type_met))[c(2,1,3)])), draw_quantiles = 0.5) + theme_bw() + theme(legend.title = element_blank()) + scale_fill_manual(values=c("#999999", "#E69F00")) + scale_colour_manual(values = c("turquoise3", "green", "red"), labels=c("5mC","5hmC","C")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/figures/20160515_combined.snv.nonsnv.cpg.M8.boxplot.nocolour.withn.pdf", width= 14/2.54, height= 10/2.54)
ggplot(snv.nonsnv.m, aes(x = factor(type_met, levels(factor(type_met))[c(2,1,3)]), y = 100*pct_met, fill = factor(ID, levels(factor(ID))[c(2,1)]))) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom") + scale_fill_manual(values=c("#999999", "#E69F00"), labels=c("SNV n=2610","non-SNV n=54479299")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()

rm(snv.nonsnv.m)



