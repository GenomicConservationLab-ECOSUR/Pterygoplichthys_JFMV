#!/bin/bash

############ Parte 1 ########## we indexed the samples

#bwa index -p 40nrPp 40nrPp_ORF.fasta.gz

############ Parte 2 ########## loops to align reads

#for i in *.fq;do
#name=`echo $i| cut -d "_" -f 1`;
#bwa mem -M -t 2 40nrPp $i > ${name}.sam;
#done

############ Parte 3 ########## loop to convert the file .sam to .bam

#for file in *.sam;
#do samtools view -@ 4 -b $file > $file.bam
#done

########### loop to sort the outputs .bam

#for file in *.bam;
#do samtools sort -@ 4 -o $file.sorted.bam $file.bam
#done

#for f in *.bam; do
#samtools sort -@ 4 -o $f.sorted.bam $f
#done


########## loop to join all files .sorted.bam into one

#for f in *.sorted.bam; do
#samtools merge -@ 4 merged.all.bam $f;
#done

############ we index the result

#samtools index merged.all.bam
#done

