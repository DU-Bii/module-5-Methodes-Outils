# Lignes de commande

## Chargement de l'environnement

```bash
$ ssh <username>@core.cluster.france-bioinformatique.fr
$ module load conda
$ source activate eba2018_variant_calling_python3
```

## Données

```bash
$ ls -lh /shared/data/projects/du_bii_2019/data/module5/seance3/fastq
$ ls -lh /shared/data/projects/du_bii_2019/data/module5/seance3/genome
$ ls -lh /shared/data/projects/du_bii_2019/data/module5/seance3/alignment_bwa
```

## Préparation de l'espace de travail

```bash
#Se déplacer dans son home
$ cd ~
#Créer un répertoire de travail
$ mkdir M5S3
#Copier les données du TP
$ cp -r /shared/data/projects/du_bii_2019/data/module5/seance3/* M5S3
#Se déplacer dans le dossier alignment_bwa
$ cd  ~/M5S3/alignment_bwa
```

## Contrôle qualité des données alignées

```bash
#Lancement de samtools
$ samtools --version # affiche la version (v.1.9)
$ samtools flagstat # affiche l'aide

$ srun samtools flagstat SRR1262731_extract.sort.bam >  SRR1262731.flagstat.txt

$ cat SRR1262731.flagstat.txt # visualisation du résultat
```

## Ajout des ReadGroups

```bash
$ samtools view # affiche l'aide
$ samtools view -H SRR1262731_extract.sort.bam  | grep “^@RG”

$ picard AddOrReplaceReadGroups --version # affiche la version (v2.18.9)
$ picard AddOrReplaceReadGroups --help # affiche l'aide

$ srun picard AddOrReplaceReadGroups I=SRR1262731_extract.sort.bam \
O=SRR1262731_extract.sort.rg.bam RGID=1 RGPL=Illumina RGPU=PU \
RGSM=SRR1262731 RGLB=LB
```

## Marquage des duplicats PCR

```bash
$ picard MarkDuplicates --help # affiche l’aide
$ srun --mem=8GB picard -Xmx8G MarkDuplicates I=SRR1262731_extract.sort.rg.bam \
O=SRR1262731_extract.sort.rg.md.bam M=SRR1262731_extract_metrics_md.txt \
VALIDATION_STRINGENCY=SILENT
$ srun samtools flagstat SRR1262731_extract.sort.rg.md.bam \
> SRR1262731_extract.md.flagstat.txt
$ cat SRR1262731_extract.md.flagstat.txt # nombre de duplicats
$ grep -A1 "LIBRARY" SRR1262731_extract_metrics_md.txt | less -S # % de pcrDup
```

## Filtres sur les alignements

```bash
# Suppression des reads non mappés et filtre sur les reads avec MAPQ < 30
$ srun samtools view -bh -F 4 -q 30 SRR1262731_extract.sort.rg.md.bam \
> SRR1262731_extract.sort.rg.md.filt.bam

$ srun samtools flagstat SRR1262731_extract.sort.rg.md.filt.bam \
> SRR1262731_extract.filt.flagstat.txt

$ cat SRR1262731_extract.filt.flagstat.txt

# Conservation des alignements dans les régions ciblées
$ bedtools --version # affiche la version (v2.27.1)
$ bedtools intersect  --help # affiche l’aide

$ srun bedtools intersect -a SRR1262731_extract.sort.rg.md.filt.bam \
-b ../additionnal_data/QTL_BT6.bed \
> SRR1262731_extract.sort.rg.md.filt.onTarget.bam

$ srun samtools index  SRR1262731_extract.sort.rg.md.filt.onTarget.bam
```

## Analyse de la couverture

```bash

# Calcul de la couverture avec samtools
$ samtools depth --help # affiche l’aide

$ srun samtools depth -b ../additionnal_data/QTL_BT6.bed \
SRR1262731_extract.sort.rg.md.filt.onTarget.bam \
> SRR1262731_extract.onTarget.depth.txt   

$ head SRR1262731_extract.onTarget.depth.txt 
```

## GATK avec sortie VCF

```bash
$ gatk HaplotypeCaller --version   # affiche la version de GATK (v 4.0.10.0)
$ gatk HaplotypeCaller  # affiche l’aide d’HaplotypeCaller
Required Arguments:
--input,-I:String         	BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
--output,-O:String        	File to which variants should be written  Required.
--reference,-R:String     	Reference sequence file  Required.

--min-base-quality-score,-mbq:Byte
                    Minimum base quality required to consider a base for calling  Default value: 10.

$ cd ../
$ mkdir -p GATK/vcf
$ cd GATK/
# Détection de variant GATK avec sortie VCF
$ srun --mem=8G gatk HaplotypeCaller \
-I ../alignment_bwa/SRR1262731_extract.sort.rg.md.filt.onTarget.bam \
-L ../additionnal_data/QTL_BT6.bed \
-O vcf/SRR1262731_extract_GATK.vcf \
-R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
--min-base-quality-score 18 \
--minimum-mapping-quality 30 \
-ERC "NONE"
```

## GATK avec sortie gVCF

```bash
# Détection de variants GATK avec sortie gVCF
$ srun --mem=8GB gatk HaplotypeCaller \
-I ../alignment_bwa/SRR1262731_extract.sort.rg.md.filt.onTarget.bam \
-L ../additionnal_data/QTL_BT6.bed -O gvcf/SRR1262731_extract_GATK.g.vcf \
-R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
--min-base-quality-score 18 \
--minimum-mapping-quality 30 \
 -ERC "GVCF"

$ ls -ltrh gvcf/

# Fusion des fichier gVCF en un seul gVCF
$ srun --mem=8GB gatk CombineGVCFs \
-L ../additionnal_data/QTL_BT6.bed \
-R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
--variant gvcf/SRR1262731_extract_GATK.g.vcf \
--variant gvcf/SRR1205992_extract_GATK.g.vcf \
--variant gvcf/SRR1205973_extract_GATK.g.vcf \
-O gvcf/pool_GATK.g.vcf

# Détection de variants simultanée sur les 3 échantillons du gVCF
$ srun --mem=8GB gatk GenotypeGVCFs -R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
--variant gvcf/pool_GATK.g.vcf \
-O gvcf/pool_GATK.vcf
```

## Samtools mpileup/varscan2 1

```bash
# Affichage de l’aide de samtools mpileup
$ samtools mpileup
Usage: samtools mpileup [options] in1.bam [in2.bam [...]]
-q, --min-MQ INT    	skip alignments with mapQ smaller than INT [0]
# L’aide de Varscan s’affiche avec le lancement de $ varscan (v2.4.3)
# Affichage de l’aide de varscan mpileup2cns
$ varscan mpileup2cns -h
USAGE: java -jar VarScan.jar mpileup2cns [pileup file] OPTIONS
mpileup file - The SAMtools mpileup file
OPTIONS:
    --min-coverage    Minimum read depth at a position to make a call [8]
    --min-reads2    Minimum supporting reads at a position to call variants [2]
    --min-avg-qual    Minimum base quality at a position to count a read [15]
```

## Samtools mpileup/varscan2 2

```bash
# Creation d’un nouveau dossier
$ cd ..
$ mkdir -p Varscan
$ cd Varscan

# Conversion du fichier d’alignement “bam” en format “mpileup”
$ srun samtools mpileup -q 30 -B -d 10000 -f ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
../alignment_bwa/SRR1262731_extract.sort.rg.md.filt.onTarget.bam \
 > SRR1262731_extract.mpileup # -A pour garder les paires anormales
 # Détection de variants avec Varscan
$ srun varscan mpileup2cns SRR1262731_extract.mpileup --output-vcf --variants --min-avg-qual 18 > SRR1262731_extract_Varscan.vcf
```

## VCF-merge

```bash
$ bgzip # v1.9
$ tabix # v1.9
$ vcftools # v0.1.16

# Renommer l’échantillon dans le VCF
$ sed -i 's|Sample1|SRR1262731.Varscan|g' SRR1262731_extract_Varscan.vcf 

# Compression et indexation du fichiers vcf
$ bgzip -c SRR1262731_extract_Varscan.vcf > SRR1262731_extract_Varscan.vcf.gz
$ tabix -p vcf SRR1262731_extract_Varscan.vcf.gz

# Merge des trois échantillons appelés avec Varscan
$ srun vcf-merge SRR1262731_extract_Varscan.vcf.gz SRR1205992_extract_Varscan.vcf.gz SRR1205973_extract_Varscan.vcf.gz > pool_Varscan.vcf
```

## SelectVariants et HardFiltering 1

```bash
# Préparation d’un nouveau répertoire de résultats
$ cd ..
$ mkdir filter_and_annot
$ cd filter_and_annot

# Extraction des SNVs dans un fichier séparé pour GATK
$ srun --mem=8GB gatk SelectVariants -R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
-V ../GATK/gvcf/pool_GATK.vcf \
-O pool_GATK.SNP.vcf \
--select-type SNP

# Extraction des SNVs dans un fichier séparé pour Varscan
$ srun --mem=8GB gatk SelectVariants -R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
-V ../Varscan/pool_Varscan.vcf \
-O pool_Varscan.SNP.vcf \
--select-type SNP
```

## SelectVariants et HardFiltering 2

```bash
# Filtrage des SNVs selon les filtres recommandés par GATK
$ srun --mem=8GB gatk VariantFiltration -R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
-V pool_GATK.SNP.vcf \
-O pool_GATK.SNP.prefilt.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name “hard_filtering_snv”


# Sélection des variants passant ce filtre
$ srun --mem=8GB gatk SelectVariants -R ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa \
-V pool_GATK.SNP.prefilt.vcf  \
-O pool_GATK.SNP.filtered.vcf  \
--exclude-filtered
```

## Intersection des résultats des variant callers

```bash
# Intersection des variants obtenus avec Varscan et avec GATK
$ vcftools # v0.1.16

# Compression et indexation des fichiers vcfs
$ srun bgzip -c pool_GATK.SNP.filtered.vcf > pool_GATK.SNP.filtered.vcf.gz
$ srun tabix -p vcf pool_GATK.SNP.filtered.vcf.gz

$ srun bgzip -c pool_Varscan.SNP.vcf > pool_Varscan.SNP.vcf.gz
$ srun tabix -p vcf pool_Varscan.SNP.vcf.gz

$ srun vcf-isec -f -n +2 pool_GATK.SNP.filtered.vcf.gz pool_Varscan.SNP.vcf.gz > GATK_varscan_inter.vcf 
```

## SnpEff

```bash
# Création de la base de données SnpEff
$ snpEff -version # affiche la version (v4.3t)

$ echo BosTaurus.genome >> snpeff.config # <genome_name>.genome 
$ mkdir -p BosTaurus
$ cp ../genome/Bos_taurus.UMD3.1.dna.toplevel.6.fa BosTaurus/sequences.fa
$ cp ../genome/Bos_taurus.UMD3.1.93_6.gtf BosTaurus/genes.gtf
$ echo -e "BosTaurus\nSnpEff4.1" > BosTaurus.db

$ srun snpEff build -c snpeff.config -gtf22 -v BosTaurus -dataDir .

# Annotation avec notre base de données
$ srun snpEff eff -c snpeff.config -dataDir . BosTaurus -s snpeff_resultat.html GATK_varscan_inter.vcf > GATK_varscan_inter.annot.vcf
```

## SnpSift

```bash
$ SnpSift filter -h  # affiche l’aide (v 4.3t)

# Garder les variants codant qui ne sont pas des synonymes :
$ srun cat GATK_varscan_inter.annot.vcf | SnpSift filter "(ANN[*].EFFECT != 'synonymous_variant') && (ANN[*].BIOTYPE = 'protein_coding')" > GATK_varscan_inter.annot.coding.nosyn.vcf

# Sélectionner notre variant d’intérêt parmi les variants hétérozygotes ayant un impact (missense) 
$ srun cat GATK_varscan_inter.annot.coding.nosyn.vcf | SnpSift filter "ANN[*].EFFECT = 'missense_variant' & isHet( GEN[2] ) & isVariant( GEN[2] ) & isRef( GEN[0] ) & isRef( GEN[1] )" > GATK_varscan_inter.annot.coding.nosyn.filtered.vcf
```

