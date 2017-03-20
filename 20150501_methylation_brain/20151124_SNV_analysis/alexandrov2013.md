
# How many raw mutations?

wc -l /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/alexandrov2013/20151119/somatic_mutation_data/Glioblastoma/Glioblastoma_raw_mutations_data.txt
# 26812

grep "Glioblastoma" /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/alexandrov2013/20151119/cancer_type_summary.txt
# Glioblastoma	 98 	 -   	 98 	 26,812 	 23,217 	 3,595 	 3,572

# https://github.com/Rdatatable/data.table/wiki
# In R:

library(data.table)
mydata<-fread("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/alexandrov2013/20151119/somatic_mutation_data/Glioblastoma/Glioblastoma_raw_mutations_data.txt")
dim(mydata)
# [1] 26812     8

length(table(mydata$V1)) # Number of Glioblastoma samples
# [1] 98

mydata[V2=="subs"]
mydata[V2=="subs", .(.N), by=.(V3,V4,V6,V7)]
table(mydata[V2=="subs", .(.N), by=.(V3,V4,V6,V7)]$N)
mydata[V2=="subs", .(.N), by=.(V3,V4,V6,V7)][N==38]
mydata[V2=="subs" & V4==124331867]

# Only C>T
mydata[V2=="subs" & V6=="C" & V7=="T"]
mydata[V2=="subs" & V6=="C" & V7=="T", .(.N), by=.(V3,V4,V6,V7)]
table(mydata[V2=="subs" & V6=="C" & V7=="T", .(.N), by=.(V3,V4,V6,V7)]$N)
mydata[V2=="subs" & V6=="C" & V7=="T", .(.N), by=.(V3,V4,V6,V7)][N==38]
mydata[V2=="subs" & V6=="C" & V7=="T" & V4==124331867]


# Only C>C
mydata[V2=="subs" & (V6=="C" & V7=="C")] # 273 examples, this must be an error


# Only C>any or any>C
mydata[V2=="subs" & (V6=="C" | V7=="C")]
mydata[V2=="subs" & (V6=="C" | V7=="C"), .(.N), by=.(V3,V4,V6,V7)]
table(mydata[V2=="subs" & (V6=="C" | V7=="C"), .(.N), by=.(V3,V4,V6,V7)]$N)
mydata[V2=="subs" & (V6=="C" | V7=="C"), .(.N), by=.(V3,V4,V6,V7)][N==38]
mydata[V2=="subs" & (V6=="C" | V7=="C") & V4==124331867]


# How many mutations match from our study match the gliblastoma mutations in Alexandrov2013?
mydata[V2=="subs"]

library(VariantAnnotation)
vcf_tvm<-readVcf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/vcf/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf", "hg19")
vcf_tvm_snv <- vcf_tvm[isSNV(vcf_tvm)]
vcf_tvb<-readVcf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/vcf/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.vcf", "hg19")
vcf_tvb_snv <- vcf_tvb[isSNV(vcf_tvb)]
vcf_mvb<-readVcf("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/vcf/CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.vcf", "hg19")
vcf_mvb_snv <- vcf_mvb[isSNV(vcf_mvb)]

tvm_snv_chrm<-as.vector(seqnames(rowRanges(vcf_tvm_snv)))
tvm_snv_start<-start(ranges(rowRanges(vcf_tvm_snv)))
tvm_snv_end<-end(ranges(rowRanges(vcf_tvm_snv)))
tvm_snv_ref<-as.vector(ref(vcf_tvm_snv))
tvm_snv_alt<-as.vector(unlist(alt(vcf_tvm_snv)))
tvm_snv<-data.table(chrm=tvm_snv_chrm, start=tvm_snv_start, end=tvm_snv_end, ref=tvm_snv_ref, alt=tvm_snv_alt)

mydata2<-mydata[V2=="subs"][,list(V3, V4, V5, V6, V7)]
setnames(mydata2, c("chrm", "start", "end", "ref", "alt"))
setkey(mydata2, chrm, ref, alt, start, end)
ov<-foverlaps(tvm_snv, mydata2, by.x=c("chrm", "ref", "alt", "start", "end"), by.y=c("chrm", "ref", "alt", "start", "end"))
match(FALSE, is.na(ov$start))
# [1] 588
ov[588,]
#    chrm ref alt    start      end  i.start    i.end
# 1:   10   C   T 73560497 73560497 73560497 73560497

# Take-home message: only 1 out of the 88 coding variants are found in the Exome sequencing results of Glioblastoma as in the Alexandrov paper.

# What do we know about that particular variant?
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
chrms<-paste0("chr", seqlevels(vcf_tvm_snv))
chrms<-replace(chrms, chrms=="chrMT", "chrM")
vcf_tvm_snv <- renameSeqlevels(vcf_tvm_snv, chrms)
intersect(seqlevels(vcf_tvm_snv), seqlevels(txdb))
seqlengths(vcf_tvm_snv)["chrM"]<-seqlengths(txdb)["chrM"]
coding <- predictCoding(vcf_tvm_snv, txdb, seqSource=Hsapiens)
unlist(lapply(split(coding$CONSEQUENCE, coding$QUERYID),unique))

coding[start(ranges(coding))==73560497,] # GENEID 64072

library(org.Hs.eg.db)

select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns=c("ENZYME")) # EC number
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="GENENAME") # Gene name
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="GO") # GO id
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="ONTOLOGY") # GO molecular function
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="PATH") # KEGG PATHWAY
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="PFAM") # PFAM id
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="PMID") # Pubmed ids
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="PROSITE") # PROSITE signatures
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="SYMBOL") # GENE SYMBOL
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="UNIGENE") # UNIGENE id
select(org.Hs.eg.db, keys=c("64072"), keytype="ENTREZID", columns="UNIPROT") # Uniprot accession id
# Q9H251 - cadherin 23



