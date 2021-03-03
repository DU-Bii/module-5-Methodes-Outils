#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB


module load star
module load subread

## Mus musculus Gencode GTF
gtf="/shared/projects/dubii2020/data/rnaseq/gtf/gencode.vM22.annotation.gtf"

## Mus musculus mm10 STAR index
index="/shared/bank/mus_musculus/mm10/star-2.7.2b/"

## Toy dataset - SRR1106775
r1="/shared/projects/dubii2020/data/rnaseq/rawdata/SRR1106775-1M_1.fastq.gz"
r2="/shared/projects/dubii2020/data/rnaseq/rawdata/SRR1106775-1M_2.fastq.gz"


##-----------------------------------
## STAR MAPPING
##-----------------------------------

starOpts="--outSAMmultNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --outSAMprimaryFlag OneBestScore --outMultimapperOrder Random --outSAMattributes All"
odir="./mapping"
prefix=$(basename $r1 | sed -e 's/.fastq.gz//')
cpus=4

mkdir -p ${odir}
STAR \
  --genomeDir ${index} \
  --sjdbGTFfile ${gtf} \
  --readFilesIn ${r1} ${r2} \
  --runThreadN ${cpus} \
  --runMode alignReads \
  --outSAMtype BAM Unsorted  \
  --readFilesCommand zcat \
  --outFileNamePrefix ${odir}/${prefix}  \
  --quantMode GeneCounts \
  --outSAMattrRGline ID:${prefix} SM:${prefix} LB:Illumina PL:Illumina  \
  ${starOpts}

##---------------------------------------
## FEATURE COUNTS
##----------------------------------------
  
odir="./counts"
bam="./mapping/SRR1106775-1M_1Aligned.out.bam"

mkdir -p ${odir}
featureCounts \
  -T ${cpus} \
  -a ${gtf} \
  -o ${odir}/counts.csv \
  -p \
  -s 0 \
  ${bam}
