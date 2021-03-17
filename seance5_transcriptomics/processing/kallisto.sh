#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB

module load kallisto

## Mus musculus Gencode GTF
gtf="/shared/projects/dubii2020/data/rnaseq/gtf/gencode.vM22.annotation.gtf"

## Mus musculus mm10 STAR index
kallistoIndex="/shared/projects/dubii2020/data/rnaseq/kallisto/gencode.vM22.transcripts_index"

## Toy dataset - SRR1106775
r1="/shared/projects/dubii2020/data/rnaseq/rawdata/SRR1106775-1M_1.fastq.gz"
r2="/shared/projects/dubii2020/data/rnaseq/rawdata/SRR1106775-1M_2.fastq.gz"


##-----------------------------------
## Kallisto
##-----------------------------------


odir="./kallisto"
cpus=4

mkdir -p ${odir}
kallisto quant \
    -i ${kallistoIndex} \
    -t ${cpus} \
    -b 100 \
    --genomebam \
    -g ${gtf} \
    -o ${odir} \
    ${r1} ${r2}
