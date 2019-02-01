## ChIP-seq : "a tutorial for DUBii_2019"
## Olivier Kirsh
## Samples from GSE53366
## this not a real bash script, more a step by step recipe for training

# requires a bash terminal, an internet connection and an nncr account

# all the uncommented line (without a #) can be runned one by one as easy as "copy paste" on a bash terminal
# don't forget srun before each command if you are on the nncr cluster

# retreive and get the data from GEODataset / SRA
# Open this link on your pc
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53366

#select SRA run selector / GEO dataset
#https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=GSE53366&go=go
#find one ZBTB33 IP sample and the corresponding input
#highlight files and click on accession list button to generate an SRR_Acc_List.txt file
SRR1176031  #input
SRR1176061  #IP TF



# Connect to nncr
ssh <login>@core.cluster.france-bioinformatique.fr

#activate conda
module load conda

#check available env
conda info -e

#Activate the env which contains most of the tools (fastcq, sratool kit, samtools, bowtie, macs, etc...)
source activate dubii2019_m5_methodes_outils_tools_r

# Check where sratoolkit prefetch function will download the sra file
vdb-config -i
# you should see somthing like
/shared/home/<yourlogin>/ncbi/public
# or
/shared/projects/du_bii_2019/<yourlogin>/ncbi/public

# check where you are
pwd
# you should see
/shared/home/<yourlogin>
# or
/shared/projects/du_bii_2019<yourlogin>

# create a project folders
mkdir -p ChIPseq/sra ChIPseq/fastq ChIPseq/bam

############
# transfert your SRR_Acc_List.txt  file to your nncr account
# go back on a terminal on your pc
scp ~/downloads/SRR_Acc_List.txt  ssh <login>@core.cluster.france-bioinformatique.fr:/shared/projects/du_bii_2019/<yourlogin>/ChIPseq
###########

# go back on your nncr terminal

# you can alternatively create a text file named SRR_Acc_List.txt on your pc or on the nncr cluster

# enter in your project folder
cd ChIPseq

nano SRR_Acc_List.txt
# add the following lines
SRR1176031
SRR1176061

# type ^x, then O , then enter

# downloaded in NCBI/public/SRA/  with sratoolkit
# https://ncbi.github.io/sra-tools/
cat SRR_Acc_List.txt | xargs -n1 prefetch

# alternatively you can type
prefetch SRR1176031
prefetch SRR1176061

# move files (from where vdb-config tells you where the files are)
mv /shared/projects/du_bii_2019<yourlogin>/ncbi/public/sra/*.sra sra/


# extract fastq with sra toolkit
# https://ncbi.github.io/sra-tools/fastq-dump.html
fastq-dump --outdir fastq/ sra/SRR1176031.sra
fastq-dump --outdir fastq/ sra/SRR1176061.sra

# renamfiles
cd fastq
mv SRR1176031.fastq input.fastq
mv SRR1176061.fastq quisuisje.fastq


# FastQC
mkdir ../qc
fastqc --outdir ../qc *.fastq
cd ..

# Trimming if necessary
# with trim_galore, trimmomatic, sickle, cutadapt.
# depends on your criteria

# Mapping on hg19 genome
# index for most mapping tools can be generated with adhoc functions or downloaded from
# https://support.illumina.com/sequencing/sequencing_software/igenome.html
# or
# ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
# hg19 indexes are stored and shared in this folder
/shared/projects/du_bii_2019/data/banks/hg19b2

#define variable for index path
b2ref=/shared/projects/du_bii_2019/data/banks/hg19b2/hg19

# mapping with Bowtie2
# http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# add srun -c 8 before the bowtie
bowtie2 -t -q -p 8 --fast --phred33 -x $b2ref fastq/input.fastq > bam/input.sam
bowtie2 -t -q -p 8 --fast --phred33 -x $b2ref fastq/quisuisje.fastq > bam/quisuisje.sam


# File conversion
cd bam
samtools view -b -h -S input.sam > input.bam
samtools view -b -h -S quisuisje.sam > quisuisje.bam

# Sort
mkdir sorted
cd sorted
samtools sort ../input.bam -o input.sorted.bam
samtools sort ../quisuisje.bam -o quisuisje.sorted.bam

# Index for IGV visualization
samtools index input.sorted.bam
samtools index quisuisje.sorted.bam

cd ..


cd ..

###### optionnal
# Filtrer pour un chromosome
samtools view -h -b input.sorted.bam "chr1" > inputChr1.bam
samtools view -h -b quisuisje.sorted.bam "chr1" > quisuisjeChr1.bam
######



# Peak calling with HOMER suite, other tools like macs2 works to
# homer require some configuration we cannot do by ourself like
# perl /shared/mfs/data/software/miniconda/envs/dubii2019_m5_methodes_outils_tools_r/share/homer-4.9.1-6/.//configureHomer.pl -install hg19
# this allow peaks annotation and motif search
# for now only hg19 is configured
# to look for other genomes, type
perl /shared/mfs/data/software/miniconda/envs/dubii2019_m5_methodes_outils_tools_r/share/homer-4.9.1-6/.//configureHomer.pl -list

# for details about homer
# http://homer.ucsd.edu/homer/
# http://homer.ucsd.edu/homer/ngs/index.html
# here we work with one replicate, homer can handle replicates, read the manual!!

# all the following command can be saved in a bash script
# and run with
srun peakcalling.sh

# create bash script
nano peakcalling.sh

# copy all the lines bellow , add the require #SBATCH options

#!/bin/sh


mkdir ana bedgraph

# make tagdirectories
# http://homer.ucsd.edu/homer/ngs/tagDir.html
# here we use default settings, but many arguments exists to enrich your outputs
# we could have writen one line per sample, but here it's a for loop


cd bam

for n in *.bam
  do 
  nom=${n%%.*}
  echo $nom
  cd ../ana
  makeTagDirectory ${nom}_tagdir/ ../bam/$nom.bam  
  cd ../bam	
done


# create Bedgraphs
for n in *.bam
  do
  nom=${n%%.*}
  echo $nom
  makeUCSCfile ../ana/${nom}_tagdir -o ../bedgraph/$nom.bedgraph -name $nom
done

# perform peak calling
cd ana
ls -d * # list and check directory names
findPeaks quisuisje_tagdir -style factor -o auto -i input_tagdir

# Transform homer output to bed file format
cd quisuisje_tagdir
pos2bed.pl peaks.txt > ../peaks.bed
cd ../

# annotation
annotatePeaks.pl   peaks.bed hg19 > annotatedpeaks.txt


# type ^x, o, enter to save your script

# make it executable
chmod 755 peakcalling.sh

# et voila
