
New runs of BS sequencing produced by Illumina - see email from Zoya Kingsbury zkingsbury@illumina.com 01/12/15 11:02

### Where is the data?

In sblab-srv001,

```bash
cd /data/sblab-data/berald01/RunData
ls -lh
```

```
total 534G
-rw-r--r-- 1 berald01 sblab 167G Dec  2 17:02 151021_HSQ9102_0379_AC6F7TANXX.tar.gz
-rw-r--r-- 1 berald01 sblab 198G Dec  2 19:04 151023_HSQ9103_0409_BC8T2UANXX.tar.gz
-rw-r--r-- 1 berald01 sblab 170G Dec  2 20:47 151027_HSQ9102_0380_BC8T60ANXX.tar.gz
```

### Pre-processing

How much space is available?

```bash
df -h
```

There are 9.1T available in /data.

Decompress, in three different tabs,

```bash
tar zxvf 151021_HSQ9102_0379_AC6F7TANXX.tar.gz
```

```bash
pigz -dc 151023_HSQ9103_0409_BC8T2UANXX.tar.gz | tar xvf -
```

```bash
pigz -dc 151027_HSQ9102_0380_BC8T60ANXX.tar.gz | tar xvf -
```


Add [Data] to the first row/first column of each SampleSheet.csv,

```bash
sed -i '1i[Data],,,,,,,,,' /data/sblab-data/berald01/RunData/illumina/upload/main/151021_HSQ9102_0379_AC6F7TANXX/SampleSheet.csv
sed -i '1i[Data],,,,,,,,,' /data/sblab-data/berald01/RunData/illumina/upload/main/151023_HSQ9103_0409_BC8T2UANXX/SampleSheet.csv
sed -i '1i[Data],,,,,,,,,' /data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/SampleSheet.csv
```

Run bcl2fastq v2.15.0.4 as installed in sblab-srv001,

```bash
cd /data/sblab-data/berald01/RunData/illumina/upload/main/151021_HSQ9102_0379_AC6F7TANXX
bcl2fastq
```

```bash
cd ../151023_HSQ9103_0409_BC8T2UANXX
bcl2fastq
```

```bash
cd ../151027_HSQ9102_0380_BC8T60ANXX
bcl2fastq
```

The three directories give an error related to unable to open positions file for a lane and tile:

```
2015-12-09 09:00:49 [7fba75146700] INFO: Opened BCL file '/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/BaseCalls/L001/C209.1/s_1_2210.bcl.gz' with 2134488 clusters
2015-12-09 09:00:49 [7fba73d44700] ERROR: Thread: 0 caught an exception first: /tmp/bcl2fastq/src/cxx/lib/data/PositionsFile.cpp(566): Throw in function void bcl2fastq::data::PositionsFile::openFile(const boost::filesystem::path&, bool, bool, bcl2fastq::data::LaneNumber, bcl2fastq::data::TileNumber, bool)
Dynamic exception type: boost::exception_detail::clone_impl<bcl2fastq::common::InputDataError>
std::exception::what: Unable to open positions file for lane #1 and tile #2210 in "/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/": aggregate_tiles_flag=NO, is_patterned_flowcell=NO

2015-12-09 09:00:49 [7fba74745700] INFO: Opened BCL file '/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/BaseCalls/L001/C206.1/s_1_2210.bcl.gz' with 2134488 clusters
2015-12-09 09:00:49 [7fba76548700] INFO: Opened BCL file '/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/BaseCalls/L001/C208.1/s_1_2210.bcl.gz' with 2134488 clusters
2015-12-09 09:00:50 [7fba5794d700] WARNING: Rethrowing a thread exception 
2015-12-09 09:00:50 [7fba5794d700] ERROR: Thread: 0 caught an exception first: /tmp/bcl2fastq/src/cxx/lib/data/PositionsFile.cpp(566): Throw in function void bcl2fastq::data::PositionsFile::openFile(const boost::filesystem::path&, bool, bool, bcl2fastq::data::LaneNumber, bcl2fastq::data::TileNumber, bool)
Dynamic exception type: boost::exception_detail::clone_impl<bcl2fastq::common::InputDataError>
std::exception::what: Unable to open positions file for lane #1 and tile #2210 in "/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/": aggregate_tiles_flag=NO, is_patterned_flowcell=NO

2015-12-09 09:00:50 [7fba7654a780] WARNING: Rethrowing a thread exception 
2015-12-09 09:00:50 [7fba7654a780] ERROR: bcl2fastq::common::Exception: 2015-Dec-09 09:00:50: No such file or directory (2): /tmp/bcl2fastq/src/cxx/lib/data/PositionsFile.cpp(566): Throw in function void bcl2fastq::data::PositionsFile::openFile(const boost::filesystem::path&, bool, bool, bcl2fastq::data::LaneNumber, bcl2fastq::data::TileNumber, bool)
Dynamic exception type: boost::exception_detail::clone_impl<bcl2fastq::common::InputDataError>
std::exception::what: Unable to open positions file for lane #1 and tile #2210 in "/data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/": aggregate_tiles_flag=NO, is_patterned_flowcell=NO
```

We have decided to run it unzip the files locally to test with the bcl2fastq v2.17.1.14 installed at martin03@nm149s012925.

First, plug in hard disk (8a14ce64-bf5d-445c-a67d-e5834570862f) and decompress within local partition, add [Data] to the first row/first column of SampleSheet.csv and run bcl2fastq,

```bash
cd /media/martin03/8a14ce64-bf5d-445c-a67d-e5834570862f/RunData
tar zxvf 151021_HSQ9102_0379_AC6F7TANXX.tar.gz -C /media/martin03/120C15AD0C158D3B/Users/Public
sed -i '1i[Data],,,,,,,,,' /media/martin03/120C15AD0C158D3B/Users/Public/illumina/upload/main/151021_HSQ9102_0379_AC6F7TANXX/SampleSheet.csv
bcl2fastq
```

It works! I stopped this local job, it seems the version installed in sblab-srv001 (v2.15.0.4) is too old, we installed v2.17.1.14 in sblab-srv001 and run everything there again. 

Back to sblab-srv001, in three different tabs,

```bash
cd /data/sblab-data/berald01/RunData/illumina/upload/main/151021_HSQ9102_0379_AC6F7TANXX
bcl2fastq
```

```bash
cd ../151023_HSQ9103_0409_BC8T2UANXX
bcl2fastq
```

```bash
cd ../151027_HSQ9102_0380_BC8T60ANXX
bcl2fastq
```

Moving fastq files to lustre,  creating folders in nas-srv001:

```bash
mkdir /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata
mkdir /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151021_HSQ9102_0379_AC6F7TANXX
mkdir /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151023_HSQ9103_0409_BC8T2UANXX
mkdir /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151027_HSQ9102_0380_BC8T60ANXX
```

Now in sblab-srv001, in three different tabs,

```bash
rsync -arvuP --remove-source-files /data/sblab-data/berald01/RunData/illumina/upload/main/151021_HSQ9102_0379_AC6F7TANXX/Data/Intensities/BaseCalls/*.fastq.gz martin03@nas-srv001:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151021_HSQ9102_0379_AC6F7TANXX
```

```bash
rsync -arvuP --remove-source-files /data/sblab-data/berald01/RunData/illumina/upload/main/151023_HSQ9103_0409_BC8T2UANXX/Data/Intensities/BaseCalls/*.fastq.gz martin03@nas-srv001:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151023_HSQ9103_0409_BC8T2UANXX
```

```bash
rsync -arvuP --remove-source-files /data/sblab-data/berald01/RunData/illumina/upload/main/151027_HSQ9102_0380_BC8T60ANXX/Data/Intensities/BaseCalls/*.fastq.gz martin03@nas-srv001:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata/151027_HSQ9102_0380_BC8T60ANXX
```


### Rename

```bash
ssh -Y uk-cri-lcst01
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata
python
```

```python
import os
import collections
prefixes= {'M8BS': 'ear042_', 'M8oxBS': 'ear043_', 'T3BS': 'ear044_', 'T3oxBS': 'ear045_', 'Undetermined': ''}
fcounter= collections.Counter()
allfiles= []
for root, dirs, files in os.walk("/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata"):
    runid= os.path.split(root)[1]
    files.sort()
    for f in files:
        if f.endswith(".fastq.gz"):
            id= f.split('_')[0]
            laneName= '_'.join(f.split('_')[1:])
            prefix= prefixes[id] + id
            newname= '.'.join([prefix, runid, laneName])
            print os.path.join(root, f), newname
            os.rename(os.path.join(root, f), newname)
            fcounter[runid]+=1
            allfiles.append(newname)


assert len(fcounter) == 3
assert len(set(allfiles)) == 80 * 3
for x in fcounter:
    assert fcounter[x] == 80


quit()
```



### Trim, align, clip

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/

ref=/lustre/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

for fq1 in fastq/20151208_BSnewdata/ear*_R1_001.fastq.gz
do
	fq2=${fq1/_R1_/_R2_}

	out=`basename $fq1 _R1_001.fastq.gz`
	out1=${out}.R1.fq.gz
	out2=${out}.R2.fq.gz

	echo $fq1, $fq2, $out, $out1, $out2

	if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
	if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi

	bsub -R "rusage[mem=20000]" -J bwameth-$out -oo bsub/$out.log "
	cutadapt -m 10 -O 1 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o fastq_trimmed/$out1 -p fastq_trimmed/$out2 $fq1 $fq2 && 
#	rsync --remove-source-files --progress -arvu $fq1 sblab-srv001:/nas/sblab_data1/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/ && 
#	rsync --remove-source-files --progress -arvu $fq2 sblab-srv001:/nas/sblab_data1/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/ &&
	bwameth.py -t 4 --reference $ref --prefix bwameth/20151208_BSnewdata/${out}.unclipped fastq_trimmed/$out1 fastq_trimmed/$out2 && 
	rm fastq_trimmed/$out1 && 
	rm fastq_trimmed/$out2 && 
	bam clipOverlap --in bwameth/20151208_BSnewdata/${out}.unclipped.bam --out bwameth/20151208_BSnewdata/${out}.bam --stats --storeOrig XC && 
	rm bwameth/20151208_BSnewdata/${out}.unclipped.bam &&
	rm bwameth/20151208_BSnewdata/${out}.unclipped.bam.bai" 
done
```


### Move .fastq files to archive

We do this after alignment (once we are confident fastq.gz are ok, we have used them and we no longer need them). Be aware there is no option to remove data from /archive directly after 48h.

```bash
ssh -Y nas-srv001
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/fastq/20151208_BSnewdata
ls -l
ls -l /archive/Groups/SBLab/fs04/martin03/repository/fastq
mv *.fastq.gz /archive/Groups/SBLab/fs04/martin03/repository/fastq/
```

### Merge and mark duplicates

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata

runBatch='hiseq201512' # id for this batch of runs to distinguish from previous batches

libids=`ls -1 ear04[2-5]*_L00[1-8].bam | cut -f1 -d'.' | sort | uniq`

ended=""
for id in $libids
do
bams=`echo $id*.S?_L00?.bam`
n=`echo $bams | wc -w`
    if [ $n != 24 ] 
    then
        echo 'WRONG'; break
    fi
    bsub $ended -J merge-${id} -R "rusage[mem=4096]" -oo ${id}.log " 
    samtools merge -f -@8 ${id}.${runBatch}.bam $bams &&
    samtools index ${id}.${runBatch}.bam &&
    java -Xmx3G -jar ~/bin/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT TMP_DIR=./ I=${id}.${runBatch}.bam O=/dev/null M=${id}.${runBatch}.markdup.txt &&
    rm $bams
    "
    ended="-w ended(merge-${id})"
done
```

Repeat all above because went over disk quota usage in picard.jar, I am going to ask for more space.

### 2016-01-19 Methylation whole genome¶

```
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/

## Split chroms in chunks to speed up and reduce memory.
ref=/lustre/sblab/berald01/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

chroms=`cut -f 1 $ref.fai`

mkdir chroms
for bam in ear042_M8BS.hiseq201512.bam ear043_M8oxBS.hiseq201512.bam ear044_T3BS.hiseq201512.bam ear045_T3oxBS.hiseq201512.bam
do
    for chr in $chroms
    do
        bdg=${bam%.bam}.$chr.bedGraph
        bsub -w "ended(merge-*)" -J $bdg -oo chroms/$bdg.log -R "rusage[mem=8192]" "bam2methylation.py -l /lustre/sblab/berald01/reference_data/hg19.allCpG.bed.gz -i $bam -r $ref -A -s ' -q 10 -F2820' -mm -mq 15 -R $chr | pigz -f > chroms/${bdg}.gz" 
    done
done
```

### Collapse CpG sites

```bash
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/chroms/

# Intersect CpG with bedgraph files, sum counts within CpGs and extract those 
# where both C+ and C- have a count
for bdg in ear04*.hiseq201512.chr*.bedGraph.gz
do
    bname=${bdg%.bedGraph.gz}
    chr=${bdg#*.chr}
    chr=chr${chr%%.*}
    echo -e " 
    tabix /lustre/sblab/berald01/reference_data/hg19.allCpG.bed.gz $chr \
    | intersectBed -a - -b $bdg -wa -wb -sorted | \
    awk -v OFS=\" \" '{print \$1, \$2, \$3, \$12, \$13, \$4}' | \
    groupBy -i - -g 6 -c 1,2,3,4,5,6 -o first,min,max,sum,sum,count | \
    awk -v OFS=\" \" '{if(\$7 > 1) print \$2, \$3, \$4, \$5, \$6}' > $bname.cpg.bdg" > $bdg.tmp.sh
    
    bsub -J cpg-${bname} -R "rusage[mem=2048]" -oo $bname.log < $bdg.tmp.sh
done

# Check each interval is exactly 2. For each file it should return nothing:
for bdg in *.hiseq201512.*.cpg.bdg
do
echo $bdg
awk '{print $3-$2}' $bdg | grep -v 2
done

## Concatenate files within libraries
for id in ear042_M8BS.hiseq201512 ear043_M8oxBS.hiseq201512 ear044_T3BS.hiseq201512 ear045_T3oxBS.hiseq201512
do
    bsub "cat ${id}.chr*.cpg.bdg | bgzip > ${id}.cpg.bedGraph.gz &&
    tabix -p bed ${id}.cpg.bedGraph.gz" 
done
rm *.chr*.cpg.bdg

## Concatenate files
mvsync *.cpg.bedGraph.gz sblab-srv001:/nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/
mvsync *.cpg.bedGraph.gz.tbi sblab-srv001:/nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/

rm -r /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/chroms/
```

### Consistency previous batch

*NB* Files from old and new batch moved to `hiseq201509-12`

```R
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/

for bdg in *.cpg.bedGraph.gz
do
tabix $bdg chr18 > /tmp/${bdg}.txt
done
tableCat.py -i /tmp/*.cpg.bedGraph.gz.txt -r '.cpg.bedGraph.gz.txt' > /tmp/chr18.tmp.bdg
rm /tmp/*.cpg.bedGraph.gz.txt

R
library(data.table)
library(ggplot2)
bdg<- fread('/tmp/chr18.tmp.bdg')
setnames(bdg, names(bdg), c('chrom', 'start', 'end', 'cnt_met', 'cnt_tot', 'ds'))
bdg[, pct_met := cnt_met/cnt_tot * 100]
met<- dcast.data.table(data= bdg, end ~ ds, value.var= 'cnt_met')
tot<- dcast.data.table(data= bdg, end ~ ds, value.var= 'cnt_tot')

ft.pval<- rep(NA, 10000)
j<- 1
for(i in seq(1, nrow(met), length.out= 10000)){
    x<- matrix(as.numeric(c(met[i, list(ear042_M8BS.hiseq201509, ear042_M8BS.hiseq201512)], tot[i, list(ear042_M8BS.hiseq201509, ear042_M8BS.hiseq201512)])), nrow= 2)
    if(any(is.na(x)) == FALSE && all(x > 20)){
        ft.pval[j]<- fisher.test(x, conf.int= FALSE)$p.value
    }
    j<- j+1
}
pdf('fisher_test_pval.pdf')
hist(ft.pval)
dev.off()

bdgct<- dcast.data.table(data= bdg, end ~ ds, value.var= 'pct_met')

ggplot(data= bdgct[seq(1, nrow(bdgct), length.out= 10000)], aes(x= ear045_T3oxBS.hiseq201509, y= ear045_T3oxBS.hiseq201512)) +
    geom_point(alpha= 0.8, colour= 'grey', size= 0.5, alpha= 0.2) +
    geom_abline(slope= 1, intercept= 0, colour= 'red', linetype= 'dashed') +
    geom_smooth(method= 'lm')
ggsave('pctmet_old_new.pdf')

smry<- bdg[, list(cnt_met= sum(cnt_met, na.rm= TRUE), cnt_tot= sum(cnt_tot, na.rm= TRUE)), by= ds]
smry[, pct_met := cnt_met/cnt_tot * 100]

#                          ds  cnt_met  cnt_tot  pct_met
#1:   ear042_M8BS.hiseq201509 17532676 24096346 72.76072
#2:   ear042_M8BS.hiseq201512 16754521 23007418 72.82226
#3: ear043_M8oxBS.hiseq201509 12906753 23481901 54.96469
#4: ear043_M8oxBS.hiseq201512 12172367 22143382 54.97068
#5:   ear044_T3BS.hiseq201509 15989177 23825394 67.10981
#6:   ear044_T3BS.hiseq201512 15395829 22918854 67.17539
#7: ear045_T3oxBS.hiseq201509 16617538 25345628 65.56373
#8: ear045_T3oxBS.hiseq201512 15197128 23185985 65.54446
```

### Combinine datasets

```bash
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/

for tbx in ear04{2..5}
do
    ls ${tbx}*.hiseq201509.cpg.bedGraph.gz ${tbx}*.hiseq201512.cpg.bedGraph.gz
    nohup zcat ${tbx}*.hiseq201509.cpg.bedGraph.gz ${tbx}*.hiseq201512.cpg.bedGraph.gz \
    | sort -k1,1 -k2,2n \
    | mergeBed -d -1 -c 4,5 -o sum,sum \
    | awk '$5 > 0' \
    | bgzip > ${tbx}.cpg.bedGraph.gz &
done

for tbx in ear04{2..5}.cpg.bedGraph.gz
do
    nohup tabix -p bed $tbx &
done
mv *.hiseq2015* hiseq201509-12/

## Rename
mv ear042.cpg.bedGraph.gz ear042_M8BS.cpg.bedGraph.gz
mv ear043.cpg.bedGraph.gz ear043_M8oxBS.cpg.bedGraph.gz
mv ear044.cpg.bedGraph.gz ear044_T3BS.cpg.bedGraph.gz
mv ear045.cpg.bedGraph.gz ear045_T3oxBS.cpg.bedGraph.gz

mv ear042.cpg.bedGraph.gz.tbi ear042_M8BS.cpg.bedGraph.gz.tbi
mv ear043.cpg.bedGraph.gz.tbi ear043_M8oxBS.cpg.bedGraph.gz.tbi
mv ear044.cpg.bedGraph.gz.tbi ear044_T3BS.cpg.bedGraph.gz.tbi
mv ear045.cpg.bedGraph.gz.tbi ear045_T3oxBS.cpg.bedGraph.gz.tbi

```

### Prepare bedgraph files for 5mC and 5hmC

Useful for deepTools. First create a bedgraph with % 5mC as value then convert bedGraph to bigWig. For 5hmC you need to subtract BS-oxBS.

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/

nohup zcat ear043_M8oxBS.cpg.bedGraph.gz | awk -v OFS=' ' '{if($5 == 0){print $0; exit} print $1, $2, $3, $4/$5 * 100}' > ear043_M8oxBS.cpg.tmp.bedGraph &
nohup zcat ear045_T3oxBS.cpg.bedGraph.gz | awk -v OFS=' ' '{if($5 == 0){print $0; exit} print $1, $2, $3, $4/$5 * 100}' > ear045_T3oxBS.cpg.tmp.bedGraph &

nohup bedGraphToBigWig ear043_M8oxBS.cpg.tmp.bedGraph /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai margin.pct_5mC.bw &
nohup bedGraphToBigWig ear045_T3oxBS.cpg.tmp.bedGraph /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai tumor.pct_5mC.bw &

## 5hmC
intersectBed -a ear042_M8BS.cpg.bedGraph.gz -b ear043_M8oxBS.cpg.bedGraph.gz -sorted -wa -wb \
| awk -v OFS=' ' '{bs=$4/$5*100; oxbs=$9/$10*100; hmc=bs - oxbs; if(hmc < 0){hmc= 0}; print $1, $2, $3, hmc}' > margin.pct_5hmC.cpg.tmp.bedGraph

intersectBed -a ear044_T3BS.cpg.bedGraph.gz -b ear045_T3oxBS.cpg.bedGraph.gz -sorted -wa -wb \
| awk -v OFS=' ' '{bs=$4/$5*100; oxbs=$9/$10*100; hmc=bs - oxbs; if(hmc < 0){hmc= 0}; print $1, $2, $3, hmc}' > tumor.pct_5hmC.cpg.tmp.bedGraph

nohup bedGraphToBigWig margin.pct_5hmC.cpg.tmp.bedGraph /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai margin.pct_5hmC.bw &
nohup bedGraphToBigWig tumor.pct_5hmC.cpg.tmp.bedGraph /data/sblab-data/common/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai tumor.pct_5hmC.bw &
```

### Plot global 5mC and 5hmC

```R
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/
R
library(data.table)
library(ggplot2)

marHmc<- fread('bigWigToBedGraph margin.pct_5hmC.bw /dev/stdout')
marHmc[, sample := 'margin']
setnames(marHmc, names(marHmc), c('chrom', 'start', 'end', 'pct_hmc', 'sample'))
tumHmc<- fread('bigWigToBedGraph tumor.pct_5hmC.bw /dev/stdout')
tumHmc[, sample := 'tumor']
setnames(tumHmc, names(tumHmc), c('chrom', 'start', 'end', 'pct_hmc', 'sample'))

hmc<- rbindlist(list(marHmc, tumHmc))
setnames(hmc, names(hmc), c('chrom', 'start', 'end', 'pct_hmc', 'sample'))
hmc[, start := NULL]

ggplot(data= NULL) +
    geom_histogram(aes(x= marHmc$pct_hmc[seq(1, nrow(marHmc), length.out= 1000000)], ..density..), fill= 'blue', alpha= 0.5) +
    geom_histogram(aes(x= tumHmc$pct_hmc[seq(1, nrow(tumHmc), length.out= 1000000)], ..density..), fill= 'red', alpha= 0.5) +
    ggtitle('5hmC in margin (blue) and tumor (red)') +
    xlim(0, 75) +
    xlab('%5hmC')
ggsave('hist_5hmC.pdf', width= 12/2.54, height= 10/2.54)
system('rsync --remove-source-files hist_5hmC.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151106_genome_methyl/')

smry<- hmc[, list(mean= mean(pct_hmc), stdev= sd(pct_hmc)), by= sample]
# sample      mean     stdev
# margin 17.494518 12.780160
#  tumor  3.422519  4.998139

## If 40% of the tumor is margin cells, how much 5hmC does the margin contributes?
17.494518 * 0.4 ~ 6.99%

## What margin contamination is needed to observe 3.4% hmC in tumor?
3.422519 / 17.494518 = 19.6%

####
marMc<- fread('bigWigToBedGraph margin.pct_5mC.bw /dev/stdout')
marMc[, sample := 'margin']
setnames(marMc, names(marMc), c('chrom', 'start', 'end', 'pct_mc', 'sample'))
tumMc<- fread('bigWigToBedGraph tumor.pct_5mC.bw /dev/stdout')
tumMc[, sample := 'tumor']
setnames(tumMc, names(tumMc), c('chrom', 'start', 'end', 'pct_mc', 'sample'))

#mc<- rbindlist(list(marMc, tumMc))
#setnames(mc, names(mc), c('chrom', 'start', 'end', 'pct_hmc', 'sample'))
#hmc[, start := NULL]

ggplot(data= NULL) +
    geom_histogram(aes(x= marMc$pct_mc[seq(1, nrow(marMc), length.out= 1000000)], ..density..), fill= 'blue', alpha= 0.5) +
    geom_histogram(aes(x= tumMc$pct_mc[seq(1, nrow(tumMc), length.out= 1000000)], ..density..), fill= 'red', alpha= 0.5) +
    ggtitle('5mC in margin (blue) and tumor (red)') +
    xlim(0, 100) +
    xlab('%5mC')
ggsave('hist_5mC.pdf', width= 12/2.54, height= 10/2.54)
system('rsync --remove-source-files hist_5mC.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151106_genome_methyl/')


```

5mC and 5hmC on aggregated counts, possibly better as it gives the correct weight to each position.

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/
for bdg in *.cpg.bedGraph.gz
do
    echo $bdg
    zcat $bdg | awk '{met+=$4; tot+=$5}END{print met, tot, NR}'
done
                               C           CT          CNT       PCT_MET  PCT_HMC
ear042_M8BS.cpg.bedGraph.gz    1335925615  1919646852  27253473  0.696
ear043_M8oxBS.cpg.bedGraph.gz  968096860   1858365119  27196199  0.521    0.175
ear044_T3BS.cpg.bedGraph.gz    1342784236  2062332798  27231244  0.651
ear045_T3oxBS.cpg.bedGraph.gz  1364404667  2149538614  27217723  0.635    0.016
```

* Expected 5hmC in tumor given 40% margin and 0 hmC in tumor:

(0.175 x 0.4) + (0 x 0.6) ~ 7% 

* What margin contamination is needed to observe 1.6% hmC in tumor?

0.016 / 0.175 ~ 9%

## Histograms of read coverage

```
cd /nas/sblab_data1/berald01/projects/20150501_methylation_brain/20151106_genome_methyl/
R
library(ggplot2)
library(data.table)

mbs<- fread('zcat ear042_M8BS.cpg.bedGraph.gz')
mox<- fread('zcat ear043_M8oxBS.cpg.bedGraph.gz')
tbs<- fread('zcat ear044_T3BS.cpg.bedGraph.gz')
tox<- fread('zcat ear045_T3oxBS.cpg.bedGraph.gz')

mbs[, sample := 'Margin BS']
mox[, sample := 'Margin oxBS']
tbs[, sample := 'Tumor BS']
tox[, sample := 'Tumor oxBS']

dat<- rbindlist(list(mbs, mox, tbs, tox))
setnames(dat, 'V5', 'tot')

gg<- ggplot(data= dat[seq(1, nrow(dat), length.out= 5000000)], aes(x= ifelse(tot > 200, 200, tot))) +
    geom_histogram(aes(y= ..density..)) +
    xlab('Depth (count of M + U) in CpG context') +
    ggtitle('Read depth') +
    facet_wrap(~sample)
ggsave('20160126_hist_cnt_tot.pdf', w= 16, h= 16, units= 'cm')
system('rsync --remove-source-files 20160126_hist_cnt_tot.pdf $mac_office:$cri_public_projects/20150501_methylation_brain/20151106_genome_methyl/')
```

## 22-03-2016 bam files moved

Top make space on `lustre` bam files from `bwameth` have been to Dario's external disk drive:

On Dario's Ubuntu desktop:

```
cd "/media/berald01/Seagate Backup Plus Drive/nas/sblab_data1/berald01/repository/bwameth/"
for bam in ear042_M8BS.hiseq201512.bam \
           ear043_M8oxBS.hiseq201512.bam \
           ear044_T3BS.hiseq201512.bam \
           ear045_T3oxBS.hiseq201512.bam 
do
rsync -arvuP --remove-source-files hpcgate.cri.camres.org:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/$bam ./
rsync -arvuP --remove-source-files hpcgate.cri.camres.org:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/$bam.bai ./
done

for bam in ear042_M8BS.hiseq201509.bam \
           ear043_M8oxBS.hiseq201509.bam \
           ear044_T3BS.hiseq201509.bam \
           ear045_T3oxBS.hiseq201509.bam
do
rsync -arvuP --remove-source-files hpcgate.cri.camres.org:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/$bam ./
rsync -arvuP --remove-source-files hpcgate.cri.camres.org:/lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/$bam.bai ./
done
```