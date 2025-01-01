#!/bin/bash

# Datos GBS (Genome by Sequencing)
# NexSeq 1x90 pb
# Enzymes used psti y mspi (clean samples)
# 103 individuos 

# Pipeline to perform quality control of reads, trimming to specific lengths, RAD locus construction, SNP calling and SNP filtering. 

###################################### Parte 1 ###########################################################################

DIRECTORIOS="Carpeta1 Lecturas2 denovo3 Poblacion4 VCFtools5"

if [ -d "$DIRECTORIOS" ]
then
   echo "Los directorios de trabajo existen"
else
    mkdir ${DIRECTORIOS} 
fi 



myPATH=`pwd`

# mv *.fq $myPATH/Lecturas2
# mv gbstrim.pl $myPATH/Carpeta1 
# mv popmap denovo1
###############################################################################################################################
# Control de Calidad
################################################################################################################################
# cd $myPATH/Carpeta1
# fastqc file.fastq


################################################################################################################################
# Trimming of adapters and length of reads
################################################################################################################################
# cd $myPATH/Carpeta1
# ls *.fastq > Lista
# for i in $(cat Lista); do name=`echo $i | cut -d "_" -f 2|cut -d "." -f 1` ; perl gbstrim.pl --enzyme1 psti --enzyme2 mspi --read R1 --minlength 50 --fastqfile $i --outputfile $myPATH/Lecturas2/$name.fastq ;done

#for infile in /home/vboxuser/Documentos/P3/ *.fq ;  do  base=$(basename --suffix=.fq $infile);  TrimmomaticSE -threads 2 \ ${infile} ${base}_trimmed.fq \
  CROP:90 MINLEN:90
 done


# fastqc *.fastq

######################################### Parte 2 ##############################################################################
################################################################################################################################ 
# STACKS::denovo (https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) 
# Run ustacks, cstacks, sstacks, tsv2bam, gstacks
# SNP calls
################################################################################################################################
# How to create a popmap?
# cd $myPATH/Lecturas2
#Crear un popmap:
# ls *.fq | sed 's/.fq//g' > popmap
# nano popmap 

# cd $myPATH/denovo3
# mv popmap $denovo3
# cd $myPATH/denovo3

 denovo_map.pl --samples $myPATH/Lecturas2/ --popmap popmap --out-path ./ -M 5 -n 5 --threads 4
# corrida population
 
# populations -P ./ -M ./popmap -r 0.8 --vcf -t 4

 # Running descriptors
# stacks-dist-extract denovo_map.log cov_per_sample > Descriptores1.txt
# stacks denovo_map --samples
# Compare optimization runs: 
# stacks-dist-extract denovo_map.log cov_per_sample > Descriptores2.txt

################################################################################################################################
# STACKS::populations
# Run to generate population-level summary statistics and export data in a variety of formats
################################################################################################################################
# cd $myPATH/Poblacion4

# populations -P $myPATH/denovo3 -O $myPATH/poblacion4 -M $myPATH/denovo3/popmap --write-random-snp --vcf 

# cd $myPATH/denovo3

# denovo_map.pl --samples $myPATH/Lecturas2/ --popmap popmap --out-path ./ -M 3 -n 3 


################################################################################################################################
# VCFtools 
# Filtrar SNP's 
################################################################################################################################

# $myPATH/VCFtools5
# vcftools --vcf $myPATH/Poblacion4/populations.snps.vcf --maf 0.01 --min-alleles  2 --max-alleles 2 --hwe 0.00001 --max-missing 0.8 --recode


