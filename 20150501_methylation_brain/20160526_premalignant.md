

The goal of this script is to test the ideas discussed with Euni and Dario about Tumour, Margin-premalignant and Margin-healthy reads.



### Where is the SNV bedGraph?

/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph

Copy it locally. In nm149s012925:

```bash
rsync martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/snv/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph /media/martin03/Seagate\ Backup\ Plus\ Drive/martin03/tmp/

```



### Where is the CpG sites bedGraph for hg19?

/lustre/sblab/martin03/reference_data/hg19.allCpG.bed.gz



### Where are the Tumour oxBS bam files? 

They were moved to Dario's local hardisk Seagate Backup Plus Drive:

https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151208_BSnewdata/20151208_BSnewdata.md#22-03-2016-bam-files-moved

I copied them to the local computer. In nm149s012925:

```bash
cp /media/martin03/Seagate\ Backup\ Plus\ Drive/nas/sblab_data1/berald01/repository/bwameth/* /media/martin03/Seagate\ Backup\ Plus\ Drive1/martin03/bwameth/

```



### Processing

We need to get the oxBS reads in tumour overlapping both SNV and CpG sites.

In nm149s012925:

```bash
cd /media/martin03/Seagate\ Backup\ Plus\ Drive/martin03


# oxBS reads aligning to first mutation
samtools view bwameth/ear045_T3oxBS.hiseq201509.bam chr1:1241100-1241100 -c # 35 hiseq201509
samtools view bwameth/ear045_T3oxBS.hiseq201512.bam chr1:1241100-1241100 -c # 39 hiseq201512


# oxBS reads aligning to SNVs
samtools view bwameth/ear045_T3oxBS.hiseq201509.bam -L tmp/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph -@ 8 -b > bwameth/ear045_T3oxBS.hiseq201509.SNV.bam
samtools view bwameth/ear045_T3oxBS.hiseq201509.SNV.bam -c # 809575
samtools view bwameth/ear045_T3oxBS.hiseq201509.SNV.bam | less
samtools view bwameth/ear045_T3oxBS.hiseq201512.bam -L tmp/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph -@ 8 -b > bwameth/ear045_T3oxBS.hiseq201512.SNV.bam
samtools view bwameth/ear045_T3oxBS.hiseq201512.SNV.bam -c # ?


# Are the bam files sorted? https://www.biostars.org/p/5256/
samtools view -H bwameth/ear045_T3oxBS.hiseq201509.bam | grep "SO" # @HD	VN:1.3	SO:coordinate
samtools view -H bwameth/ear045_T3oxBS.hiseq201512.bam | grep "SO" # @HD	VN:1.3	SO:coordinate
# both are sorted


# merge hiseq201509 and hiseq201512 and create index for it
samtools merge bwameth/ear045_T3oxBS.bam bwameth/ear045_T3oxBS.hiseq201509.bam bwameth/ear045_T3oxBS.hiseq201512.bam
samtools index bwameth/ear045_T3oxBS.bam
samtools view bwameth/ear045_T3oxBS.bam chr1:1241100-1241100 -c # hopefully hiseq201509 + hiseq201512


# visualise the output bam with igv

intersectBed -abam bwameth/ear045_T3oxBS.hiseq201512.bam -b tmp/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.bedGraph | samtools view -c

```



