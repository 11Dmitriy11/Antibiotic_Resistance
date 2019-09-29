#!/bin/bash
path1=./Reference_genome 
if [[ ! -e $path1 ]]
then
    mkdir Reference_genome 
fi
cd ./Reference_genome
fna=./GCA_000005845.2_ASM584v2_genomic.fna
if [[ ! -f $fna ]]
then
    echo "Downloading reference genom GCA_000005845.2_ASM584v2.fna GCA_000005845.2_ASM584v2"
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz
    gunzip GCA_000005845.2_ASM584v2_genomic.fna.gz 
fi
gff=./GCA_000005845.2_ASM584v2_genomic.gff
if [[ ! -f $gff ]]
then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.gff.gz
    gunzip GCA_000005845.2_ASM584v2_genomic.gff.gz 
fi
cd ../
assembly_genome=./assembly_genome
if [[ ! -e $assembly_genome ]]
then
    mkdir assembly_genome
fi
cd ./assembly_genome

forward_as=./amp_res_1.fastq
if [[ ! -f $forward_as ]]
then
    echo "Downloading amp_res assembly"
    wget http://public.dobzhanskycenter.ru/mrayko/amp_res_1.fastq.zip
    unzip amp_res_1.fastq.zip
fi
reversed_as=./amp_res_2.fastq
if [[ ! -f $reversed_as ]]
then
    wget http://public.dobzhanskycenter.ru/mrayko/amp_res_2.fastq.zip
    unzip amp_res_2.fastq.zip
fi
date +%F >> ../laboratory_book.txt
echo -n "Number of reads forvard: " >> ../laboratory_book.txt
echo "$(cat amp_res_1.fastq | wc -l)/4" | bc >> ../laboratory_book.txt
echo -n "Number of reads reversed: " >> ../laboratory_book.txt
echo "$(cat amp_res_2.fastq | wc -l)/4" | bc >> ../laboratory_book.txt

cd ../

fastqc=./fastqc_results
if [[ ! -e $fastqc ]]
then
    mkdir fastqc_results
fi

ampfastq1=./fastqc_results/amp_res_1_fastqc.html
ampfastq2=./fastqc_results/amp_res_2_fastqc.html
if [[ ! -f $ampfastq1 ]]
then
    fastqc -o ./fastqc_results ./assembly_genome/amp_res_1.fastq 
fi

if [[ ! -f $ampfastq2 ]]
then
    fastqc -o ./fastqc_results ./assembly_genome/amp_res_2.fastq 
fi

trimmomatic=./assembly_genome/trimmomatic
if [[ ! -e $trimmomatic ]]
then
    mkdir ./assembly_genome/trimmomatic
fi


java -jar $4 PE -phred33 ./assembly_genome/amp_res_1.fastq ./assembly_genome/amp_res_2.fastq ${trimmomatic}/_1P$1_$2.fq ${trimmomatic}/_1U$1_$2.fq ${trimmomatic}/_2P$1_$2.fq ${trimmomatic}/_2U$1_$2.fq LEADING:$1 TRAILING:$2 SLIDINGWINDOW:10:20 MINLEN:20
echo  "Number of reads after filte: " >> ./laboratory_book.txt
echo -n "Number of reads forvard: " >> ./laboratory_book.txt
echo "$(cat $trimmomatic/_1P$1_$2.fq | wc -l)/4" | bc >> ./laboratory_book.txt
echo -n "Number of reads reversed: " >> ./laboratory_book.txt
echo "$(cat $trimmomatic/_2P$1_$2.fq | wc -l)/4" | bc >> ./laboratory_book.txt

#_1P= ${trimmomatic}/_1P$1_$2.fq
#_2P= ${trimmomatic}/_2P$1_$2.fq


trimmomatic_f=./fastqc_results/trimmomatic
if [[ ! -e $trimmomatic_f ]]
then
    mkdir ./fastqc_results/trimmomatic
fi

trimfastq1=./fastqc_results/trimmomatic/_1P$1_$2_fastqc.html
trimfastq2=./fastqc_results/trimmomatic/_2P$1_$2_fastqc.html
if [[ ! -f $trimfastq1 ]]
then
    fastqc -o ./fastqc_results/trimmomatic ./assembly_genome/trimmomatic/_1P$1_$2.fq
fi

if [[ ! -f $trimfastq2 ]]
then
    fastqc -o ./fastqc_results/trimmomatic ./assembly_genome/trimmomatic/_2P$1_$2.fq
fi

bwa index ./Reference_genome/GCA_000005845.2_ASM584v2_genomic.fna

sam=./alignment_$1_$2.sam
if [[ ! -f $sam ]]
then
    bwa mem ./Reference_genome/GCA_000005845.2_ASM584v2_genomic.fna ./assembly_genome/trimmomatic/_1P$1_$2.fq ./assembly_genome/trimmomatic/_2P$1_$2.fq > alignment_$1_$2.sam
fi

bam=./alignment_$1_$2.bam
if [[ ! -f $bam ]]
then
    samtools view -S -b alignment_$1_$2.sam > alignment_$1_$2.bam
echo cat "$(samtools flagstat alignment_$1_$2.bam)" | sed -n 5p >> ./laboratory_book.txt
fi

bam_sort=./alignment_$1_$2_sorted.bam
if [[ ! -f $bam_sort ]]
then
    samtools sort alignment_$1_$2.bam > alignment_$1_$2_sorted.bam
    samtools index alignment_$1_$2_sorted.bam
fi

mpi=./my.mpileup
if [[ ! -f $mpi ]]
then
    samtools mpileup -f ./Reference_genome/GCA_000005845.2_ASM584v2_genomic.fna ./alignment_$1_$2_sorted.bam > my.mpileup
fi

vcf=./varscan_results_$3.vcf
if [[ ! -f $vcf ]]
then
    java -jar $5 mpileup2snp my.mpileup --min_var_freq $3 --variants --output-vcf 1 > varscan_results_$3.vcf
fi

