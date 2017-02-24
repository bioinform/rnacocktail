#!/bin/bash

mkdir small_example
cd small_example

echo "Download reference genome (chromosome 21) FASTA file"
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz

echo "Unzip reference genome (chromosome 21) FASTA file"
gunzip Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz

echo "Download reference annotation GTF file"
wget ftp://ftp.ensembl.org/pub/release-75//gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

echo "Unzip reference annotation GTF file"
gunzip Homo_sapiens.GRCh37.75.gtf.gz

echo "Restrict GTF to chromosome 21"
less Homo_sapiens.GRCh37.75.gtf |awk '{if ($1==21) print}' > Homo_sapiens.GRCh37.75.chromosome.21.gtf


echo "Download HISAT2 binaries"
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip

echo "Unzip HISAT2 binaries"
unzip hisat2-2.0.4-Linux_x86_64.zip

echo "Index genome with HISAT2"
./hisat2-2.0.4/hisat2-build Homo_sapiens.GRCh37.75.dna.chromosome.21.fa Homo_sapiens.GRCh37.75.dna.chromosome.21.HISAT2

echo "Test alignment step using HISAT2"
run_rnacocktail.py align --align_idx Homo_sapiens.GRCh37.75.dna.chromosome.21.HISAT2 --outdir out --workdir work --ref_gtf Homo_sapiens.GRCh37.75.chromosome.21.gtf --1 ../seq_1.fq.gz  --2 ../seq_2.fq.gz --hisat2 hisat2-2.0.4/hisat2 --hisat2_sps hisat2-2.0.4/hisat2_extract_splice_sites.py  --samtools samtools --sample A

echo "Test alignment step using HISAT2"
run_rnacocktail.py align --align_idx Homo_sapiens.GRCh37.75.dna.chromosome.21.HISAT2 --outdir out --workdir work --ref_gtf Homo_sapiens.GRCh37.75.chromosome.21.gtf --1 ../seq_1.fq.gz  --2 ../seq_2.fq.gz --hisat2 hisat2-2.0.4/hisat2 --hisat2_sps hisat2-2.0.4/hisat2_extract_splice_sites.py  --samtools samtools --sample A

echo "Download StringTie binaries"
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.4.Linux_x86_64.tar.gz

echo "Untar StringTie binaries"
tar -xzvf stringtie-1.2.4.Linux_x86_64.tar.gz

echo "Test reconstruction step using StringTie"
run_rnacocktail.py reconstruct --alignment_bam work/hisat2/A/alignments.sorted.bam --outdir out --workdir work --ref_gtf Homo_sapiens.GRCh37.75.chromosome.21.gtf --stringtie stringtie-1.2.4.Linux_x86_64/stringtie --sample A
