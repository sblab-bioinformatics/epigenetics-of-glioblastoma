
Gathering data to submit to ArrayExpress using [Annotare](https://www.ebi.ac.uk/fg/annotare/).

There are four samples we have used in the study (more were obtained from the patient but not used `/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/samples/20151204brain_sequencing_project_samples.xlsx`):

* T3 = F3 (tumour) : DNA-seq, (ox)BS-seq, RNA-seq
* M8 (margin) : DNA-seq, (ox)BS-seq
* M3 (margin) : RNA-seq
* blood: DNA-seq


## Raw

BS/oxBS-seq fastq files. In nas-srv001:

```bash
ls /archive/Groups/SBLab/fs04/martin03/repository/fastq/ear042* | wc -l # M8 BS-seq, 48 = 3 (runs) x 8 (lanes) x 2 (pair end)
ls /archive/Groups/SBLab/fs04/martin03/repository/fastq/ear043* | wc -l # M8 oxBS-seq, 48 = 3 (runs) x 8 (lanes) x 2 (pair end)
ls /archive/Groups/SBLab/fs04/martin03/repository/fastq/ear044* | wc -l # T3 BS-seq, 48 = 3 (runs) x 8 (lanes) x 2 (pair end)
ls /archive/Groups/SBLab/fs04/martin03/repository/fastq/ear045* | wc -l # T3 oxBS-seq, 48 = 3 (runs) x 8 (lanes) x 2 (pair end)
# Total: 192 files

```

RNA-seq fastq files. In nas-srv001:

```bash
ls /archive/Groups/SBLab/fs04/berald01/repository/fastq/ear047* | wc -l # F3, 4 = 2 (lanes) x 2 (pair end)
ls /archive/Groups/SBLab/fs04/berald01/repository/fastq/ear049* | wc -l # M3, 4 = 2 (lanes) x 2 (pair end)
# Total: 8 files

```

DNA-seq bam files. There are no fastq files because the Illumina pipeline did not generate them.

First of all, copy the bam files from uk-cri-lcst01 to the archive (delete them from uk-cri-lcst01 once the glioblastoma paper is more advanced). In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs05/martin03/repository/bam
cp /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_A01.bam* . # blood, .bam and .bam.bai
cp /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_C01.bam* . # margin (M8), .bam and .bam.bai
cp /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly/LP2000729-DNA_E01.bam* . # tumour (T3), .bam and .bam.bai
# Total: 3 files (.bam only)

```



## Processed

BS/oxBS-seq bedGraph files (might need to change them to .txt extension). In nas-srv001:

```bash
/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear042_M8BS.cpg.bedGraph.gz # M8 BS-seq, 1
/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear043_M8oxBS.cpg.bedGraph.gz # M8 oxBS-seq, 1
/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear044_T3BS.cpg.bedGraph.gz # T3 BS-seq, 1
/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear045_T3oxBS.cpg.bedGraph.gz # T3 oxBS-seq, 1
# Total: 4 files

```

RNA-seq transcript tables (might need to change them to .txt extension). In nas-srv001:

```bash
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq/ear047_F3.tx_quant.tsv # F3, 1
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq/ear049_M3.tx_quant.tsv # M3, 1
# Total: 2 files
```

DNA-seq vcf files (might need to change them to .txt extension). In nas-srv001:

```bash
# SNV
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.vcf.gz # T3 vs. blood, 1
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.vcf.gz # M8 vs. blood, 1
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/20151029_tumour_vs_margin/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01/SomaticVariations/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf.gz # T3 vs. M8, 1

# CNV
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz # T3 vs. blood, 1
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz # M8 vs. blood, 1
/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/cnv/20160211/LP2000729-DNA_C01_LP2000729-DNA_E01_G1_P1.somatic.SV.vcf # T3 vs. M8, 1 (contains both CNVs and SVs). CNVs correspond to lines containing Canvas

# SV
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz # T3 vs. blood, 1
/data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz # M8 vs. blood, 1
# T3 vs. M8, same as CNV above (LP2000729-DNA_C01_LP2000729-DNA_E01_G1_P1.somatic.SV.vcf) but SVs correspond lines not containing Canvas

# Total : 8 files

```

Overall we have a total of 217 raw and processed files to submit to ArrayExpress.



## Calculate Nominal Length and Nominal SDev

DNA-seq bam files. In nas-srv001:

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/Assembly
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=LP2000729-DNA_A01.bam O=LP2000729-DNA_A01_insert_size_metrics.txt H=LP2000729-DNA_A01_insert_size_histogram.pdf STOP_AFTER=100000000 # blood
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=LP2000729-DNA_C01.bam O=LP2000729-DNA_C01_insert_size_metrics.txt H=LP2000729-DNA_C01_insert_size_histogram.pdf STOP_AFTER=100000000 # margin (M8)
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=LP2000729-DNA_E01.bam O=LP2000729-DNA_E01_insert_size_metrics.txt H=LP2000729-DNA_E01_insert_size_histogram.pdf STOP_AFTER=100000000 # tumour (T3)

```

BS/oxBS-seq bam files. The original files are in `martin03@nm149s012925:/media/martin03/Seagate Backup Plus Drive/martin03/bwameth`. However there is a copy of the `hiseq201509.chr18` data in `/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth`. We will use these to calculate the insert size. In nas-srv001:

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=ear042_M8BS.hiseq201509.chr18.bam O=ear042_M8BS.hiseq201509.chr18_insert_size_metrics.txt H=ear042_M8BS.hiseq201509.chr18_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # M8 BS-seq
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=ear043_M8oxBS.hiseq201509.chr18.bam O=ear043_M8oxBS.hiseq201509.chr18_insert_size_metrics.txt H=ear043_M8oxBS.hiseq201509.chr18_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # M8 oxBS-seq
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=ear044_T3BS.hiseq201509.chr18.bam O=ear044_T3BS.hiseq201509.chr18_insert_size_metrics.txt H=ear044_T3BS.hiseq201509.chr18_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # T3 BS-seq
java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=ear045_T3oxBS.hiseq201509.chr18.bam O=ear045_T3oxBS.hiseq201509.chr18_insert_size_metrics.txt H=ear045_T3oxBS.hiseq201509.chr18_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # T3 oxBS-seq

```

We only have RNA-seq fastq files, we need to get the bam files in order to calculate the insert size. Using kallisto in nas-srv001:

Copy files.

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/rnaseq/data
mkdir fastq
cd fastq

cp /archive/Groups/SBLab/fs04/berald01/repository/fastq/ear047* .
cp /archive/Groups/SBLab/fs04/berald01/repository/fastq/ear049* .
```

Align files and produce bams following [this](https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160303_rnaseq/20160303_rnaseq.md#transcript-quantification-with-kallisto) and [this](https://pachterlab.github.io/kallisto/pseudobam.html).

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/rnaseq/data/fastq
mkdir ../bam

bsub -R "rusage[mem=8192]" -q interactive -Is bash

fq1=ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_R1_001.fastq.gz
fq2=ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_R2_001.fastq.gz

kallisto quant --bias --pseudobam -i /lustre/sblab/berald01/reference_data/kallisto/Homo_sapiens.GRCh37.cdna.all.idx -o ../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001.tx <(zcat $fq1) <(zcat $fq2) | samtools view -Sb - > ../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001.bam

java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001.bam O=../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_insert_size_metrics.txt H=../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # ear047_F3

cat ../bam/ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_insert_size_metrics.txt

fq1=ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_R1_001.fastq.gz
fq2=ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_R2_001.fastq.gz

kallisto quant --bias --pseudobam -i /lustre/sblab/berald01/reference_data/kallisto/Homo_sapiens.GRCh37.cdna.all.idx -o ../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001.tx <(zcat $fq1) <(zcat $fq2) | samtools view -Sb - > ../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001.bam

java -Xmx3G -jar ~/sw/picard/picard-tools-1.138/picard.jar CollectInsertSizeMetrics I=../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001.bam O=../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_insert_size_metrics.txt H=../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_insert_size_histogram.pdf VALIDATION_STRINGENCY=SILENT # ear049_M3

cat ../bam/ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_insert_size_metrics.txt

```


## Calculating sizes of file:

### Raw

BS/oxBS-seq fastq files. In nas-srv001:

```bash
du -ch /archive/Groups/SBLab/fs04/martin03/repository/fastq/ear04[2-5]* | grep "total" # 667G
```

RNA-seq fastq files. In nas-srv001:

```bash
du -ch /archive/Groups/SBLab/fs04/berald01/repository/fastq/ear04[79]* | grep "total" # 11G
```

DNA-seq bam files. In nas-srv001:

```bash
du -ch /archive/Groups/SBLab/fs05/martin03/repository/bam/LP2000729-DNA_[ACE]01.bam* | grep "total" # 485G
```

### Processed

BS/oxBS-seq bedGraph files. In nas-srv001:

```bash
ls -lh /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/ear04[2-5]_*.cpg.bedGraph.gz | grep -v "hiseq" # ~ 1.6G
```

RNA-seq transcript tables. In nas-srv001:

```bash
du -ch /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq/ear04[79]_*.tx_quant.tsv | grep "total" # 14M
```

DNA-seq vcf files in nas-srv001 should also be pretty light.



## Calculating md5 hashes for files

### Raw

BS/oxBS-seq fastq files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs04/martin03/repository/fastq/
md5sum ear04[2-5]*
```

154c98571b0f3eaa25778a70e6d325a4  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz
73438f3edc9c12e4308a548af6fecfff  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R2_001.fastq.gz
908f1e74540cb0049186e66f8e4cf2d3  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L002_R1_001.fastq.gz
041323c19ea71830f40371444eb33ee8  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L002_R2_001.fastq.gz
90ad0dea53e3b7eec82f12d4cb084801  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L003_R1_001.fastq.gz
49ec0707ad88107bcbfb46d14f24e8f5  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L003_R2_001.fastq.gz
51108d3608ee8feef9efb69d59fd701b  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L004_R1_001.fastq.gz
ecf7f531964aeb2d7c6d09ff7f9e9f9d  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L004_R2_001.fastq.gz
4f256f2c8e67083d68afa62cd47dd151  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L005_R1_001.fastq.gz
bce6d044657dc236e30f508d6eb619e5  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L005_R2_001.fastq.gz
771aa83883b4a3af051a46fc87b71b3c  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L006_R1_001.fastq.gz
efc75c3863d9bfa82335b15e71570d35  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L006_R2_001.fastq.gz
d63490d7b30488acfd21522c672cb773  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L007_R1_001.fastq.gz
c6eafa94d6875347ccf52d1f62e04db5  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L007_R2_001.fastq.gz
818e28fa4917a0f804814fef0e4d5cf2  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L008_R1_001.fastq.gz
6d637917d1d1a1fd5e9d2c33a45787ad  ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L008_R2_001.fastq.gz
61b256b511c8c5529461d072f1d48a04  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L001_R1_001.fastq.gz
da92851076f3372aa2f68c717c2cd2e6  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L001_R2_001.fastq.gz
65c5fb612106e9ca443883909050e63e  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L002_R1_001.fastq.gz
466152b3aee45413120535530f5d8195  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L002_R2_001.fastq.gz
6b666873afc66ad75ff7beb739d79743  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L003_R1_001.fastq.gz
1de9b5460c8c642f54aac694617a3e96  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L003_R2_001.fastq.gz
769fafe5d3809382fba7a3fee8817f3a  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L004_R1_001.fastq.gz
1c93a176216a45104cc87aa7080fa424  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L004_R2_001.fastq.gz
652942707c31430c2d1bb347663a0792  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L005_R1_001.fastq.gz
7d3e4d16844355f8dae08e8ab8c6771d  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L005_R2_001.fastq.gz
28803fe7a5a162a946ae69328bd6be3a  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L006_R1_001.fastq.gz
569bf9550b56fb86a94dc435e5a0ea56  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L006_R2_001.fastq.gz
3a18aa8c9fbd3917cee09815cb6f4de7  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L007_R1_001.fastq.gz
83ccd067764503bbb82ec6e54f48a2aa  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L007_R2_001.fastq.gz
65e8b435cc9a52e8f832498e9801e3cf  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L008_R1_001.fastq.gz
f21c4428a448e918f629cd117135aaf6  ear042_M8BS.151023_HSQ9103_0409_BC8T2UANXX.S3_L008_R2_001.fastq.gz
bf43b12141a639f01259f5ef2ca365c8  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L001_R1_001.fastq.gz
61a32c2531e95d2613209905e69f22c8  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L001_R2_001.fastq.gz
f98722eca79fa00a916d8453f5e67b6f  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L002_R1_001.fastq.gz
973759842284db3372d5af56c3009b7d  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L002_R2_001.fastq.gz
47cb27218f54b928e18cf12fff73a775  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L003_R1_001.fastq.gz
32717bc69ba8b18de7d3bb8e24092fcd  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L003_R2_001.fastq.gz
d0037f46f08de5c2e4a73c4c5b4f9b32  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L004_R1_001.fastq.gz
c681cff878ee9a4d0e193f620862cb6b  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L004_R2_001.fastq.gz
9d3bbc697e46be7bdae18f41aa25acf1  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L005_R1_001.fastq.gz
e2496a8e9d43459e6f2804cff524a539  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L005_R2_001.fastq.gz
851a1056381ed7ace8a5c3d93f145660  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L006_R1_001.fastq.gz
757b0ea5d81eff791a402039ccb9fe85  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L006_R2_001.fastq.gz
637c24b101a6600ee5850811c51f148e  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L007_R1_001.fastq.gz
3c258a177c32ef996fcfea49f524b2f9  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L007_R2_001.fastq.gz
e3cdaf11c119e0a8d211bb4755b96e91  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L008_R1_001.fastq.gz
fa84711c78680a56e3760542ed9e6198  ear042_M8BS.151027_HSQ9102_0380_BC8T60ANXX.S3_L008_R2_001.fastq.gz
9aa44cd3b065fa2567db211f05fa48db  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L001_R1_001.fastq.gz
ecefe3c03e953c4df9fda66b500884bb  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L001_R2_001.fastq.gz
c5f718c79e414eb32178e0ac43bf63a1  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L002_R1_001.fastq.gz
1a7eb3011b8ce86be96d3589d45ff5a4  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L002_R2_001.fastq.gz
c497c00e01a7741aa0c31885c3ecb058  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L003_R1_001.fastq.gz
2e4d0bcdb4f62b5fbc8813abd93ab9ea  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L003_R2_001.fastq.gz
2b5c9768f792e8c41693f19f1f04f2f6  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R1_001.fastq.gz
fa2b9bb742ffed7b7e9acefde867ffdf  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
a9e3c937d84b5eb5a9ded67891667d47  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L005_R1_001.fastq.gz
ea046c839f838b58737d615a69fe7571  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L005_R2_001.fastq.gz
1545168c7a47a0d2759eac8ddea32059  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L006_R1_001.fastq.gz
e0d50a23c17383654b785eceb4b798bf  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L006_R2_001.fastq.gz
6b9f792e9a34a1abd4250eb4d27c069b  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L007_R1_001.fastq.gz
4359864a823207a31db9339c11b0fbfd  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L007_R2_001.fastq.gz
c95a1dee25fc3cedc73ecd34a64b3a9c  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L008_R1_001.fastq.gz
d4cff99759c060ab24704a9c86c9b7fb  ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L008_R2_001.fastq.gz
11f860df96e02b3745d5d5e53011bfb1  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L001_R1_001.fastq.gz
c11f5b8fa0cf72b2e7e4d65a5f81e4ef  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L001_R2_001.fastq.gz
93f5772b94f3ee07a4bd59a8f6251cd3  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L002_R1_001.fastq.gz
6612b7decd48b9d7d8d3c7586d87947b  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L002_R2_001.fastq.gz
4fbfa6386d2645379d2d43c16491b39b  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L003_R1_001.fastq.gz
ead42a2b046153d81d53382b897f09ce  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L003_R2_001.fastq.gz
993f668d1ae82bb0880095424e894c2a  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L004_R1_001.fastq.gz
77cb266e69f954b808ac683f00d91545  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L004_R2_001.fastq.gz
435ed66d921030688bd1ebdf74fe038e  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L005_R1_001.fastq.gz
cbc1fb36b23bdc2f6c6bd20522e43639  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L005_R2_001.fastq.gz
8193aedb65dcdfd5bd9f295f5cd86fef  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L006_R1_001.fastq.gz
552f532d3c94371b9db26f12f26cb78f  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L006_R2_001.fastq.gz
3131535786f70ec5e09a42e5313e24a0  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L007_R1_001.fastq.gz
757d147f3c83fb5538c0aa3539c9e72c  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L007_R2_001.fastq.gz
08805a1b5c872e7bbf7ccf4790de2cb4  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L008_R1_001.fastq.gz
b8753dcdacf17532aabc54979c47ace7  ear043_M8oxBS.151023_HSQ9103_0409_BC8T2UANXX.S4_L008_R2_001.fastq.gz
630ff6e681838e007f296221abbe586d  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L001_R1_001.fastq.gz
f253bc840b1ad328bd09de39e759edb0  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L001_R2_001.fastq.gz
e58efd9d7d7841a059fd0dfff3aca56e  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L002_R1_001.fastq.gz
635c99fdf17be72d4914808e66c089f4  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L002_R2_001.fastq.gz
a9baee071d47334d54f42edc42cf31fe  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L003_R1_001.fastq.gz
9fe78a02adee7abfce850906d9ee96b7  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L003_R2_001.fastq.gz
8f044ab9054be0829c22022d4d1e9075  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L004_R1_001.fastq.gz
83aa1ddae20e00c2e7d418bfa0902c0c  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L004_R2_001.fastq.gz
5cb953d31f0177c5afc69326ea7fe5c4  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L005_R1_001.fastq.gz
1a36524cac60cbc813ced66b78628296  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L005_R2_001.fastq.gz
989dbfd8cb8ca37383373ba436ee4a29  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L006_R1_001.fastq.gz
e30766dde4225b8ac3e5eb1a39186242  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L006_R2_001.fastq.gz
866367105d9102877ba8b31b702e91fa  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L007_R1_001.fastq.gz
abd4b83177ad047e3b2d58894f91b989  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L007_R2_001.fastq.gz
0754ec14cf7f77e76364333a2c920953  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L008_R1_001.fastq.gz
49336503b3babe16a3c67915a480a911  ear043_M8oxBS.151027_HSQ9102_0380_BC8T60ANXX.S4_L008_R2_001.fastq.gz
c8e519e557fa1415b4f669a501cc2dd2  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L001_R1_001.fastq.gz
20c4d41ec586bc8f52ba28d133ed7e86  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L001_R2_001.fastq.gz
de9cf9e5842a03d6f2e074bd26ce8c98  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L002_R1_001.fastq.gz
e477a2b5e3abef488278c7021496605c  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L002_R2_001.fastq.gz
28f8bf32ce9278b47aa140efc3186f76  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L003_R1_001.fastq.gz
f611ed8748257c8a44519d25845d1ced  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L003_R2_001.fastq.gz
a184b7e31ceee6dcc26e56ab895c32fa  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L004_R1_001.fastq.gz
c2347ed60039c973c3a4f6a3075adccb  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L004_R2_001.fastq.gz
677b9be25b95daf5697c8d5f9d3475b6  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L005_R1_001.fastq.gz
f98314097d856116b1d226aa830e3dd5  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L005_R2_001.fastq.gz
1a029221f26f7d0d2288c2153c646f8b  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L006_R1_001.fastq.gz
a74f96dd5d00f7aeb8492789d8214765  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L006_R2_001.fastq.gz
71b35babd06dd8f7a61687ed3f1ae72a  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L007_R1_001.fastq.gz
c8eacda869f4e845b364527a40767944  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L007_R2_001.fastq.gz
17d19051f26eab67b5ac3f2682cfc952  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L008_R1_001.fastq.gz
b2054f6f6dab311768f5915008d25b85  ear044_T3BS.151021_HSQ9102_0379_AC6F7TANXX.S1_L008_R2_001.fastq.gz
5ea42873b2efc418e1d9b4f20fd84314  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L001_R1_001.fastq.gz
dc0da06dadcf7fd28a73ebe35b4db77c  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L001_R2_001.fastq.gz
4b467c39974e42361e4e86f6198e7309  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L002_R1_001.fastq.gz
ce995ddd9a62b0aa41beeea207680126  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L002_R2_001.fastq.gz
f156b7c47558bfcfb76e915d0704db57  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L003_R1_001.fastq.gz
7e89ce1407180dae3df8b49da9699b50  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L003_R2_001.fastq.gz
ea69d63e870603772fec498f139669b3  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L004_R1_001.fastq.gz
ee05818047e796fc894c81fe61bc553e  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L004_R2_001.fastq.gz
93effa0ed734581ab631070e7fcd061d  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L005_R1_001.fastq.gz
2e6a242ec9f4881070cf7c1900afadcf  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L005_R2_001.fastq.gz
a34b019991050071c936a843cbf213fd  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L006_R1_001.fastq.gz
a5bea9fd46b5099c60d0870c17e97de2  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L006_R2_001.fastq.gz
2a5b1f93a2946c12a0981830c42ff7a9  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L007_R1_001.fastq.gz
5098e4d1ecc1bb81e4ff30995a056494  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L007_R2_001.fastq.gz
c24bac3f006ec2a34b0b926f2218275c  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L008_R1_001.fastq.gz
6fa25b15cd18578681e7290e2a4ea4b5  ear044_T3BS.151023_HSQ9103_0409_BC8T2UANXX.S1_L008_R2_001.fastq.gz
215bb2b1cd395b3d43ef1947c5f66736  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L001_R1_001.fastq.gz
5a8f0690e1afa8cea0bced0b6a5fcab4  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L001_R2_001.fastq.gz
d2f9c627be5c4f0263b1aca98e52e0fb  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L002_R1_001.fastq.gz
2389d67c07e84ab99c5395d8ce27b65a  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L002_R2_001.fastq.gz
6e3efc35d1d91c98abba16c092e7e8f0  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L003_R1_001.fastq.gz
3fecd6ea4ff3db4540518537c51de423  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L003_R2_001.fastq.gz
7c165bab5c951b416398fd2c1107a62e  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L004_R1_001.fastq.gz
59dadff23b4ebae8fd7ab5d5aacfe3ca  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L004_R2_001.fastq.gz
79f0b3c5267f9300ce0e8872760cdf2d  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L005_R1_001.fastq.gz
7ec8ea5695384398c9150136a95365d5  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L005_R2_001.fastq.gz
22d5a4c850d2f532a8f0daef073fcfaa  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L006_R1_001.fastq.gz
615c74d6a0e0595259254b798cc56de6  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L006_R2_001.fastq.gz
113e21dfa846c92747f0a36df31d6a79  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L007_R1_001.fastq.gz
69f63c4d67a989dfcc18008ecf6148a5  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L007_R2_001.fastq.gz
8bcdd44f9678917bddf308bd4671bf0b  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L008_R1_001.fastq.gz
35036b12311cf6602e6805e5ce771a02  ear044_T3BS.151027_HSQ9102_0380_BC8T60ANXX.S1_L008_R2_001.fastq.gz
7709d45be176164212c65ab31ab4f132  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L001_R1_001.fastq.gz
2bdff8686f1cd68ef4261dc3b9d35c9a  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L001_R2_001.fastq.gz
004a5b7cc56d875603922ed104948bdd  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L002_R1_001.fastq.gz
bb26a77debbe1e8d4c6f6671518e5845  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L002_R2_001.fastq.gz
7db6e7aef12e8d8603aae994d69f518b  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L003_R1_001.fastq.gz
b56767d507e88d8399a16127bf8040c6  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L003_R2_001.fastq.gz
7388a23cfff58cae323443037d1cea80  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L004_R1_001.fastq.gz
cdc1c7a4ccd7290bcae097d67dd11424  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L004_R2_001.fastq.gz
c0f3ed8be9eefe8fa6ad66fed3a70710  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L005_R1_001.fastq.gz
2d92fbb7eb628d2cbe2e4034aed07c43  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L005_R2_001.fastq.gz
3924b03ff4910387cdfbaba90ff60e8b  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L006_R1_001.fastq.gz
d359cd7eee6b723d708d22f8c06a1b6a  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L006_R2_001.fastq.gz
aaf8eefb5c68ee3bcb381887cc9f4933  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L007_R1_001.fastq.gz
bd6585f4811a08624cae1a1c384aed9c  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L007_R2_001.fastq.gz
9837cf4684e9f23820dbf52047ccd47f  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L008_R1_001.fastq.gz
745bf0c8e775197ab08e2d7feab4acf6  ear045_T3oxBS.151021_HSQ9102_0379_AC6F7TANXX.S2_L008_R2_001.fastq.gz
77d7e976ef1789e8553411bddde290c0  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L001_R1_001.fastq.gz
e44c0f4554ef2c0e3e13e2e75a211823  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L001_R2_001.fastq.gz
59d4a5ca7b81e2d895fed48750b6fb3c  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L002_R1_001.fastq.gz
5a30bce07ae4374de69d1ee2cad2a45b  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L002_R2_001.fastq.gz
877868f6d88e5adc5d6f6866171e4039  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L003_R1_001.fastq.gz
6676fa521a68d9b03ee74c3a62e15b40  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L003_R2_001.fastq.gz
1d0062f270c32e46e5f9f1c6b0f34265  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L004_R1_001.fastq.gz
4c1a5686aa9bf24618991b479d0d53aa  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L004_R2_001.fastq.gz
0db3d305e93d9cd67606cec5d7bd4b70  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L005_R1_001.fastq.gz
0e534310664e8c81f4495dc96f571eee  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L005_R2_001.fastq.gz
aad29fb36f762efeb323af26a5ac06df  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L006_R1_001.fastq.gz
55615bfca46450e34f4f1ff8163db413  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L006_R2_001.fastq.gz
a1983c486aac11160a6695247bec5930  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L007_R1_001.fastq.gz
c9479c1307595bce99fec5bbc712f14f  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L007_R2_001.fastq.gz
dc27968d0c6219d8e9f5c470d09f75c9  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L008_R1_001.fastq.gz
c9afa840863e9c8cc480f4ed4845059d  ear045_T3oxBS.151023_HSQ9103_0409_BC8T2UANXX.S2_L008_R2_001.fastq.gz
a494d4477cfef02d4c86328d9b3b5e1c  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L001_R1_001.fastq.gz
6c184d6c05872a13f2577b22489c03d9  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L001_R2_001.fastq.gz
7a0ffc734bbcda18601aa62c2fa82d08  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L002_R1_001.fastq.gz
b37690658dd885f5e40b2e199bb2b367  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L002_R2_001.fastq.gz
293854a68612a9c49833e975c6842e17  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L003_R1_001.fastq.gz
05162e00c1f5cdc373ae3445f9c51aa8  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L003_R2_001.fastq.gz
c4c13eb4a392dfe9ac4d1ee133e120c4  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L004_R1_001.fastq.gz
f34fd7f25bf306f32e09d0814138fb43  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L004_R2_001.fastq.gz
d6f24e55d71c47fa18cd67006ad7c2a9  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L005_R1_001.fastq.gz
de7f5212e224abbca93a0f911d47010e  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L005_R2_001.fastq.gz
375de50fa71ae1c453318d6527910300  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L006_R1_001.fastq.gz
b806a384aa96e69f4831c41a48efe4a0  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L006_R2_001.fastq.gz
bdd372a871cf96cd2334f45b92269f5e  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L007_R1_001.fastq.gz
c938594e0d9c92a7f63053a136f83665  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L007_R2_001.fastq.gz
2f0edb4d55766678e039bc98c4dc0494  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L008_R1_001.fastq.gz
67cbca9749209809d76da5231a99d254  ear045_T3oxBS.151027_HSQ9102_0380_BC8T60ANXX.S2_L008_R2_001.fastq.gz

RNA-seq fastq files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs04/berald01/repository/fastq/
md5sum ear04[79]*
```

10ca666ea6a2ad890fcc33aecc3a297b  ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_R1_001.fastq.gz
cc98d1ddd8f9fff87113b29567b1370d  ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L001_R2_001.fastq.gz
9c40dd622ab3792d17dfce009e2380c6  ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L002_R1_001.fastq.gz
49e1faee32e0df13c189dd5aa4350647  ear047_F3.CT17574-RNA-access-CRUK-F3_S3_L002_R2_001.fastq.gz
385cbce717f4bd12995a42897d385638  ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_R1_001.fastq.gz
0ef57aac092791607ad4c68819c7cd55  ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L001_R2_001.fastq.gz
2a9c1387ecbe437d86b5e9b319b6b9cc  ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L002_R1_001.fastq.gz
8b443dfe563f20259058b69a32d6f26f  ear049_M3.CT17694-RNA-access-CRUK-M3_S4_L002_R2_001.fastq.gz

DNA-seq bam files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs05/martin03/repository/bam/
md5sum LP2000729-DNA_[ACE]01.bam*
```

ae8605f38f095dbda42042cbe9a0f215  LP2000729-DNA_A01.bam
6ab97d8cedc6504dfd8b1fd928bff07a  LP2000729-DNA_A01.bam.bai
ab307a3234efe1d75f21914b180cd802  LP2000729-DNA_C01.bam
40d1cd36f17d82995289dd00749935b2  LP2000729-DNA_C01.bam.bai
b4711bbf6fc680c3bde6f5345fca6086  LP2000729-DNA_E01.bam
e9df751c9c6b9cebc186fae1604ebae9  LP2000729-DNA_E01.bam.bai


### Processed

BS/oxBS-seq bedGraph files. In nas-srv001:

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/
md5sum ear042_M8BS.cpg.bedGraph.gz
md5sum ear043_M8oxBS.cpg.bedGraph.gz
md5sum ear044_T3BS.cpg.bedGraph.gz
md5sum ear045_T3oxBS.cpg.bedGraph.gz
```

df59dbf7a0fc9833b747e168f3e96b8a  ear042_M8BS.cpg.bedGraph.gz
8f33c0aa56c572767819358fa58badab  ear043_M8oxBS.cpg.bedGraph.gz
0614e9f3731a56abb9c254d68ba578f5  ear044_T3BS.cpg.bedGraph.gz
03db77b5ebe2199236274c2914cc9808  ear045_T3oxBS.cpg.bedGraph.gz

RNA-seq transcript tables. In nas-srv001:

```bash
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq
md5sum ear04[79]_*.tx_quant.tsv
```

8f03088f00dee02e4794f288146f43b4  ear047_F3.tx_quant.tsv
d1785fbe9b7e5f4dfd4b41320500ddc1  ear049_M3.tx_quant.tsv

DNA-seq vcf files. In nas-srv001:

```bash
# SNV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.vcf.gz
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.vcf.gz
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/20151029_tumour_vs_margin/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01/SomaticVariations/
md5sum CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf.gz

# CNV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/cnv/20160211/
md5sum LP2000729-DNA_C01_LP2000729-DNA_E01_G1_P1.somatic.SV.vcf

# SV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
md5sum CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz
```

8fcfa2d2dbe60220e76aa73ed7280a25  CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.vcf.gz
a70d6e5d8005935792ab677ac7d1704d  CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.vcf.gz
30eeec9e6626952cbe9af7527ffd74ba  CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf.gz
e56badd0beaae1ec7d683eda6163df08  CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
a508f8ed3c9caafa823e76fae7f39e55  CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
5784f27d5119e1815f16c0c8b245b85e  LP2000729-DNA_C01_LP2000729-DNA_E01_G1_P1.somatic.SV.vcf
28c0309c099e6caa1b144124c78c501d  CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz
8921074ea8e5a6363022e4b2db44c557  CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz



## Uploading files to ArrayExpress using FTP

### Raw

BS/oxBS-seq fastq files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs04/martin03/repository/fastq/

ftp ftp-private.ebi.ac.uk/it4dkwf9-35v6eih5we0y9
#ftp: ftp-private.ebi.ac.uk/it4dkwf9-35v6eih5we0y9: Name or service not known

ftp ftp://aexpress:aexpress1@ftp-private.ebi.ac.uk/it4dkwf9-35v6eih5we0y9/
#ftp: ftp://aexpress:aexpress1@ftp-private.ebi.ac.uk/it4dkwf9-35v6eih5we0y9/: Name or service not known

ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
#Connected to ftp-private.ebi.ac.uk (193.62.194.179)
#220-Private FTP
#220
#Name (ftp-private.ebi.ac.uk:martin03): aexpress
#331 Please specify the password.
#Password:
#230 Login successful.
#Remote system type is UNIX.
#Using binary mode to transfer files.

ftp> cd it4dkwf9-35v6eih5we0y9
#250 Directory successfully changed.

ftp> pwd
#257 "/it4dkwf9-35v6eih5we0y9"

ftp> status
#Connected to ftp-private.ebi.ac.uk.
#No proxy connection.
#Mode: stream; Type: binary; Form: non-print; Structure: file
#Verbose: on; Bell: off; Prompting: on; Globbing: on
#Store unique: off; Receive unique: off
#Case: off; CR stripping: on
#Ntrans: off
#Nmap: off
#Hash mark printing: off; Use of PORT cmds: on
#Tick counter printing: off

ftp> put ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz
#local: ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz remote: ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz
#27 Entering Passive Mode (193,62,194,179,176,204).
#150 Ok to send data.
#226 Transfer complete.
#3405403245 bytes sent in 450 secs (7570.51 Kbytes/sec)

ftp> ls
#227 Entering Passive Mode (193,62,194,179,169,175).
#150 Here comes the directory listing.
#-rw-rw----    1 ftp      ftp      3405403245 Sep 26 12:57 ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz
#226 Directory send OK.

ftp> delete ear042_M8BS.151021_HSQ9102_0379_AC6F7TANXX.S3_L001_R1_001.fastq.gz
#227 Entering Passive Mode (193,62,194,179,165,105).
#150 Here comes the directory listing.
#226 Directory send OK.

ftp> prompt
#Interactive mode off.

ftp> mput ear042* # Terminal tab 1
ftp> mput ear043* # Terminal tab 2
ftp> mput ear044* # Terminal tab 3
ftp> mput ear045* # Terminal tab 4

ftp> quit
```

>Dear Sergio,
>
>The MD5 checksum ensures that both copies of the file are exactly the same. An
>MD5 check error indicates that there may be some data loss during transit. In
>the FTP folder, there seems to be only file
>(ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz) which
>is of 3565112128 bytes. As a first step, I suggest verifying the file size of
>your local copy. You can also try uploading the file again if file sizes
>doesn't match. Please let us know if this doesn't work.
>
>Best regards,
>Awais
>
>Sergio Martinez Cuesta <Sergio.MartinezCuesta@cruk.cam.ac.uk> wrote:
>
>Dear Amy,
>
> After uploading the files using FTP I get an "MD5 check error" in one of
> the fastq.gz files (see screenshot attached). So I am not able to assign
> this file to the samples.
>
> I ran md5sum again, I got the same hash as expected, however I don't
> understand why annotare raises an error.
>
> fa2b9bb742ffed7b7e9acefde867ffdf 
> ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
>
> Can you fix this please?
>
> Thanks,
>
> Sergio

There is indeed a difference in sizes 3762274880 (local) vs. 3565112128 (remote). So I will try uploading it again. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs04/martin03/repository/fastq/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
ls
#227 Entering Passive Mode (193,62,194,179,158,190).
#150 Here comes the directory listing.
#-rw-rw----    1 ftp      ftp      3565112128 Sep 26 15:52 ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
#226 Directory send OK.
delete ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
put ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
ls
#227 Entering Passive Mode (193,62,194,179,188,202).
#150 Here comes the directory listing.
#-rw-rw----    1 ftp      ftp      3762274880 Sep 28 17:36 ear043_M8oxBS.151021_HSQ9102_0379_AC6F7TANXX.S4_L004_R2_001.fastq.gz
#226 Directory send OK.
quit
```


RNA-seq fastq files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs04/berald01/repository/fastq/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
mput ear04[79]*
quit
```

DNA-seq bam files. In nas-srv001:

```bash
cd /archive/Groups/SBLab/fs05/martin03/repository/bam/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
mput LP2000729-DNA_A01.bam* # Terminal tab 1
mput LP2000729-DNA_C01.bam* # Terminal tab 2
mput LP2000729-DNA_E01.bam* # Terminal tab 3
quit
```

The ftp upload of `.bam` files never really completed (they seem to have been uploaded fine anyway according to the size of the files - might be corrupted though) so I quitted all three separate terminals with `Cntr+Z` and upload the corresponding `.bam.bai` separately:

```bash
cd /archive/Groups/SBLab/fs05/martin03/repository/bam/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put LP2000729-DNA_A01.bam.bai
put LP2000729-DNA_C01.bam.bai
put LP2000729-DNA_E01.bam.bai
quit
```


### Processed

BS/oxBS-seq bedGraph files. In nas-srv001:

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
mput ear042_M8BS.cpg.bedGraph.gz ear043_M8oxBS.cpg.bedGraph.gz ear044_T3BS.cpg.bedGraph.gz ear045_T3oxBS.cpg.bedGraph.gz
quit
```

RNA-seq transcript tables. In nas-srv001:

```bash
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20160303_rnaseq
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
mput ear04[79]_*.tx_quant.tsv
quit
```

DNA-seq vcf files. In nas-srv001:

```bash
# SNV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.vcf.gz
quit

cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.vcf.gz
quit

cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/20151029_tumour_vs_margin/CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_E01_NormalLP2000729-DNA_C01.somatic.vcf.gz
quit

# CNV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
quit

cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.CNV.vcf.gz
quit

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/cnv/20160211/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put LP2000729-DNA_C01_LP2000729-DNA_E01_G1_P1.somatic.SV.vcf
quit

# SV
cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_E01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_E01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz
quit

cd /data/sblab/group_folders/berald01/projects/20150501_methylation_brain/20150914_archive_snp/disk1/CancerLP2000729-DNA_C01_NormalLP200/SomaticVariations/
ftp ftp-private.ebi.ac.uk # user: aexpress, password: aexpress1
cd it4dkwf9-35v6eih5we0y9
prompt
put CancerLP2000729-DNA_C01_NormalLP2000729-DNA_A01.somatic.SV.vcf.gz
quit

```
