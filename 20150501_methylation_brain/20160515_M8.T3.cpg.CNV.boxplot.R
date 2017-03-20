library(data.table)
library(ggplot2)

ear042_M8BS.cpg.CNV <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/data/cnv/20160515_ear042_M8BS.cpg.CNV.bedGraph")
setnames(ear042_M8BS.cpg.CNV, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "cnv_type", "cnv_number"))

ear043_M8oxBS.cpg.CNV <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/data/cnv/20160515_ear043_M8oxBS.cpg.CNV.bedGraph")
setnames(ear043_M8oxBS.cpg.CNV, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "cnv_type", "cnv_number"))

ear044_T3BS.cpg.CNV <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/data/cnv/20160515_ear044_T3BS.cpg.CNV.bedGraph")
setnames(ear044_T3BS.cpg.CNV, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "cnv_type", "cnv_number"))

ear045_T3oxBS.cpg.CNV <- fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/data/cnv/20160515_ear045_T3oxBS.cpg.CNV.bedGraph")
setnames(ear045_T3oxBS.cpg.CNV, c("chr", "start", "end", "cnt_met", "cnt_tot", "strand", "cnv_type", "cnv_number"))

setkey(ear042_M8BS.cpg.CNV, chr, start, end, strand, cnv_type, cnv_number)
setkey(ear043_M8oxBS.cpg.CNV, chr, start, end, strand, cnv_type, cnv_number)
setkey(ear044_T3BS.cpg.CNV, chr, start, end, strand, cnv_type, cnv_number)
setkey(ear045_T3oxBS.cpg.CNV, chr, start, end, strand, cnv_type, cnv_number)

M8.cpg <- merge(ear042_M8BS.cpg.CNV, ear043_M8oxBS.cpg.CNV, suffixes = c(".BS", ".oxBS"))
T3.cpg <- merge(ear044_T3BS.cpg.CNV, ear045_T3oxBS.cpg.CNV, suffixes = c(".BS", ".oxBS"))

rm(ear042_M8BS.cpg.CNV, ear043_M8oxBS.cpg.CNV, ear044_T3BS.cpg.CNV, ear045_T3oxBS.cpg.CNV)

setkey(M8.cpg, chr, start, end, strand, cnv_type, cnv_number)
setkey(T3.cpg, chr, start, end, strand, cnv_type, cnv_number)
M8.T3.cpg <- merge(M8.cpg, T3.cpg, suffixes = c(".M8", ".T3"))
rm(M8.cpg, T3.cpg)

M8.T3.cpg[cnv_type == ".", cnv_type:="REF"]
M8.T3.cpg[, pct_5mC.M8 := 100*cnt_met.oxBS.M8/cnt_tot.oxBS.M8]
M8.T3.cpg[, pct_5hmC.M8 := 100*(cnt_met.BS.M8/cnt_tot.BS.M8 - cnt_met.oxBS.M8/cnt_tot.oxBS.M8)]
M8.T3.cpg[, pct_C.M8 := 100 - 100*cnt_met.BS.M8/cnt_tot.BS.M8]
M8.T3.cpg[, c("cnt_met.BS.M8", "cnt_tot.BS.M8", "cnt_met.oxBS.M8", "cnt_tot.oxBS.M8"):= NULL]
M8.T3.cpg[, pct_5mC.T3 := 100*cnt_met.oxBS.T3/cnt_tot.oxBS.T3]
M8.T3.cpg[, pct_5hmC.T3 := 100*(cnt_met.BS.T3/cnt_tot.BS.T3 - cnt_met.oxBS.T3/cnt_tot.oxBS.T3)]
M8.T3.cpg[, pct_C.T3 := 100 - 100*cnt_met.BS.T3/cnt_tot.BS.T3]
M8.T3.cpg[, c("cnt_met.BS.T3", "cnt_tot.BS.T3", "cnt_met.oxBS.T3", "cnt_tot.oxBS.T3"):= NULL]

M8.T3.cpg.m <- reshape2::melt(M8.T3.cpg, id.vars = c("chr", "start", "end", "strand", "cnv_type", "cnv_number"), variable.name = "type_met", value.name = "pct_met")
rm(M8.T3.cpg)

#M8.T3.cpg.m[type_met %like% "M8"] # select M8 only

#levels(factor(M8.T3.cpg.m$cnv_type))

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/figures/20160515_M8.cpg.CNV.pdf", width= 14/2.54, height= 10/2.54)
ggplot(M8.T3.cpg.m[type_met %like% "M8"], aes(x = factor(type_met), y = pct_met, fill = factor(cnv_type, levels(factor(cnv_type))[c(1,4,3,2)]))) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(legend.title = element_blank()) + scale_fill_manual(values=c("darkcyan", "white", "darkorange3", "dodgerblue3")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()

pdf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/summary/figures/20160515_T3.cpg.CNV.pdf", width= 14/2.54, height= 10/2.54)
ggplot(M8.T3.cpg.m[type_met %like% "T3"], aes(x = factor(type_met), y = pct_met, fill = factor(cnv_type, levels(factor(cnv_type))[c(1,4,3,2)]))) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(legend.title = element_blank()) + scale_fill_manual(values=c("darkcyan", "white", "darkorange3", "dodgerblue3")) + xlab("") + ylab("% Modification") + scale_x_discrete(labels = c("5mC","5hmC","C")) + coord_cartesian(ylim = c(0, 100))
dev.off()


