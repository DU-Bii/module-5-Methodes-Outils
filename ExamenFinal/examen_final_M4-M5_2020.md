# DU-Bii 2020 :  Examen final Modules 4 et 5 

Nous allons vous faire faire un analyse simple, de données de reséquençage d'un génome bactérien.
Les données sont issues de cet article :  "Complete Genome Sequences of 13 Bacillus subtilis Soil Isolates for Studying Secondary Metabolite Diversity"  (doi:10.1128/MRA.01406-19)

## Objectif :
Nous alons vous demander de faire une première analyse de ces données, et de nous la rendre sous la forme d'un rapport qui trace l'ensemble des étapes suivies. 
Ce rapport devra être mis à nôtre disposition dans un dépôt public GitHub. Les analyses devront pouvoir être rejouées sur le cluster de l'IFB.

Données d'entrées :
- Identifiant du run : SRR10390685
- Génome de référence : NC_000964
    - Gff https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz
    - Fasta https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz

## Consignes :

Détailler les différentes étapes dans un rapport HTML ou PDF (généré via un Rmd ou autre).

- Téléchargement des données depuis les banques publiques
- Contrôle qualité des données brutes (reads)
    - La qualité des bases vous paraît-elle satisfaisante ? Pourquoi ?
    - Quelle est la profondeur de séquençage (calculée par rapport à la taille du génome de référence) ?
- Nettoyage des reads
    - Quel pourcentage de reads sont filtrés et pourquoi ?
- Alignement des reads contre le génome de reférence
    - Quel est le % de reads pairés alignés ?
- Extraire dans un fichier BAM les reads chevauchant à au moins 50% le gène trmNF

## Informations devant figurer dans le rapport

- Présentation (par exemple à l'aide de la commande `tree`) de l'organisation du repretoire du projet
- Justification des paramètres utilisés
- Analyse succinte des résultats obtenus après chaque outil lancé (figures, tableaux ou texte)

> N'oubliez pas les informations nécessaires à la reproductibilité des analyses !!

<!--
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz
gzip -d GCF_000009045.1_ASM904v1_genomic.fna.gz
gzip -d GCF_000009045.1_ASM904v1_genomic.gff.gz
module load sra-tools
srun --cpus-per-task 8 fasterq-dump -S -p SRR10390685 --outdir . --threads 8
pigz SRR10390685_1.fastq
pigz SRR10390685_2.fastq 
module load fastqc
mkdir QC
srun --cpus-per-task 8 fastqc SRR10390685_2.fastq.gz -o QC/ -t 8
module load fastp
srun --cpus-per-task 8 fastp --in1 SRR10390685_1.fastq.gz --in2 SRR10390685_2.fastq.gz -l 100 --out1 SRR10390685_1.cleaned.fastq.gz --out2 SRR10390685_2.cleaned.fastq.gz --unpaired1 SRR10390685_singletons.fastq.gz --unpaired2 SRR10390685_singletons.fastq.gz -w 1 -h fastp.html -t 8
module load samtools
samtools faidx GCF_000009045.1_ASM904v1_genomic.fna 
more GCF_000009045.1_ASM904v1_genomic.fna.fai 
module load bwa
srun bwa index GCF_000009045.1_ASM904v1_genomic.fna
srun --cpus-per-task=4 bwa mem GCF_000009045.1_ASM904v1_genomic.fna SRR10390685_1.cleaned.fastq.gz SRR10390685_2.cleaned.fastq.gz -t 3 | samtools view -hbS - > SRR10390685.bam
samtools flagstat SRR10390685.bam 
samtools sort SRR10390685.bam -o SRR10390685_sorted.bam


grep trmNF GCF_000009045.1_ASM904v1_genomic.gff | awk '$3=="gene"' > trmNF.gff3
module load bedtools
srun samtools index SRR10390685_sorted.bam 
bedtools intersect -a SRR10390685_sorted.bam -b trmNF.gff3 -f 0.5 > SRR10390685_on_trmNF.bam
samtools view -c SRR10390685_on_trmNF.bam

-->