
The methylation files need to be regenerated so that cnt_met and cnt_tot are calculated for C and G in an individual basis.

Possibilities:
- CpG
- non-CpG # this may have to be addressed in a later stage because bam2methylation.py can't really handle this as it stands.

Adapting from here:

https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20151208_BSnewdata/20151208_BSnewdata.md#2016-01-19-methylation-whole-genome






### Methylation CpGs in libraries hiseq201509 and hiseq201512

```bash


cp /lustre/sblab/berald01/reference_data/hg19.allCpG.bed.gz /lustre/sblab/martin03/reference_data




# hiseq201512

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/

# Split chroms in chunks to speed up and reduce memory.
ref=/lustre/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
#samtools faidx $ref
chroms=`cut -f 1 $ref.fai`

mkdir chroms_allmeth_cpg
for bam in ear042_M8BS.hiseq201512.bam ear043_M8oxBS.hiseq201512.bam ear044_T3BS.hiseq201512.bam ear045_T3oxBS.hiseq201512.bam
do
    for chr in $chroms
    do
        bdg=${bam%.bam}.$chr.bedGraph
        bsub -J $bdg -oo chroms_allmeth_cpg/$bdg.log -R "rusage[mem=16384]" "
		bam2methylation.py -l /lustre/sblab/martin03/reference_data/hg19.allCpG.bed.gz -i $bam -r $ref -A -s ' -q 10 -F2820' -mm -mq 15 -R $chr | pigz -f > chroms_allmeth_cpg/${bdg}.gz
		" 
    done
done


# In order to also consider non-CpG sites I tried without the -l option for 8GB and 16 GB but they crashed
# I am trying just in chr1 with 32GB as a test:
mkdir chroms_allmeth
bam=ear042_M8BS.hiseq201512.bam
chr=chr1
bdg=${bam%.bam}.$chr.bedGraph
bsub -J $bdg -oo chroms_allmeth/$bdg.log -R "rusage[mem=32768]" "
bam2methylation.py -i $bam -r $ref -A -s ' -q 10 -F2820' -mm -mq 15 -R $chr | pigz -f > chroms_allmeth/${bdg}.gz
"

# This is why hg19.allCpG.bed.gz was used in bam2methylation.py above



# hiseq201509

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth

mkdir chroms_allmeth_cpg
for bam in ear042_M8BS.hiseq201509.bam ear043_M8oxBS.hiseq201509.bam ear044_T3BS.hiseq201509.bam ear045_T3oxBS.hiseq201509.bam
do
    for chr in $chroms
    do
        bdg=${bam%.bam}.$chr.bedGraph
        bsub -J $bdg -oo chroms_allmeth_cpg/$bdg.log -R "rusage[mem=16384]" "
		bam2methylation.py -l /lustre/sblab/martin03/reference_data/hg19.allCpG.bed.gz -i $bam -r $ref -A -s ' -q 10 -F2820' -mm -mq 15 -R $chr | pigz -f > chroms_allmeth_cpg/${bdg}.gz
		" 
    done
done


# Check completion

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/chroms_allmeth_cpg
grep "Successfully completed." ear04*hiseq201512.chr*.bedGraph.log | wc -l # must sum up to 100
cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/chroms_allmeth_cpg
grep "Successfully completed." ear04*hiseq201509.chr*.bedGraph.log | wc -l # must also sum up to 100

```




### Concatenate chromosome files within each library hiseq201512 and hiseq201509

```bash

# hiseq201512

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/20151208_BSnewdata/chroms_allmeth_cpg

for id in ear042_M8BS.hiseq201512 ear043_M8oxBS.hiseq201512 ear044_T3BS.hiseq201512 ear045_T3oxBS.hiseq201512
do
    bsub "zcat ${id}.chr*.bedGraph.gz | bgzip > ${id}.cpg.bedGraph.gz &&
    tabix -p bed ${id}.cpg.bedGraph.gz" 
done

zcat ear042_M8BS.hiseq201512.cpg.bedGraph.gz | wc -l # 26235249
zcat ear043_M8oxBS.hiseq201512.cpg.bedGraph.gz | wc -l # 54486364
zcat ear044_T3BS.hiseq201512.cpg.bedGraph.gz | wc -l # 54564169
zcat ear045_T3oxBS.hiseq201512.cpg.bedGraph.gz | wc -l # 54534404

# Repeat ear042_M8BS.hiseq201512 because something might have gone wrong:

id=ear042_M8BS.hiseq201512
bsub "zcat ${id}.chr*.bedGraph.gz | bgzip > ${id}.cpg.bedGraph.gz &&
tabix -p bed ${id}.cpg.bedGraph.gz" 

zcat ear042_M8BS.hiseq201512.cpg.bedGraph.gz | wc -l # 54602370 - better

#rm *.chr*.bedGraph.gz
mkdir /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210
cp *.cpg.bedGraph.gz* /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/


# hiseq201509

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/bwameth/chroms_allmeth_cpg

for id in ear042_M8BS.hiseq201509 ear043_M8oxBS.hiseq201509 ear044_T3BS.hiseq201509 ear045_T3oxBS.hiseq201509
do
    bsub "zcat ${id}.chr*.bedGraph.gz | bgzip > ${id}.cpg.bedGraph.gz &&
    tabix -p bed ${id}.cpg.bedGraph.gz" 
done

zcat ear042_M8BS.hiseq201509.cpg.bedGraph.gz | wc -l # 54630318
zcat ear043_M8oxBS.hiseq201509.cpg.bedGraph.gz | wc -l # 54524016
zcat ear044_T3BS.hiseq201509.cpg.bedGraph.gz | wc -l # 54591435
zcat ear045_T3oxBS.hiseq201509.cpg.bedGraph.gz | wc -l # 54581654

# rm *.chr*.bedGraph.gz
cp *.cpg.bedGraph.gz* /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/
```




### Combine libraries hiseq201512 and hiseq201509

```bash

cd /lustre/sblab/martin03/repository/20150921_BrainMethylomeRoadMap/methylation_cpg/20160210/

for tbx in ear042_M8BS ear043_M8oxBS ear044_T3BS ear045_T3oxBS
do
    ls ${tbx}.hiseq201509.cpg.bedGraph.gz ${tbx}.hiseq201512.cpg.bedGraph.gz
    bsub -J $tbx -oo ${tbx}.cpg.bedGraph.log "zcat ${tbx}.hiseq201509.cpg.bedGraph.gz ${tbx}.hiseq201512.cpg.bedGraph.gz \
    | sort -k1,1 -k2,2n \
    | mergeBed -d -1 -c 5,6,7 -o sum,sum,distinct \
    | bgzip > ${tbx}.cpg.bedGraph.gz"
done




```















