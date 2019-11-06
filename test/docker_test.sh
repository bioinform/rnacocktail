#!/bin/bash
set -ex

# mkdir example
# chmod -R 777 example/
cd example

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Initialization"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Download reference genome (chromosome 21) FASTA file"
# wget ftp://ftp.ensembl.org/pub/release-90//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

# echo "Unzip reference genome (chromosome 21) FASTA file"
# gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

# echo "Download reference annotation GTF file"
# wget ftp://ftp.ensembl.org/pub/release-90//gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz

# echo "Unzip reference annotation GTF file"
# gunzip Homo_sapiens.GRCh38.90.gtf.gz

# echo "Restrict GTF to chromosome 21"
# less Homo_sapiens.GRCh38.90.gtf |awk '{if ($1==21) print}' > Homo_sapiens.GRCh38.90.chromosome.21.gtf

# echo "Prepare reference transcriptome FASTA file"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 gffread \
# 				/work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				-g /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
# 				-w /work_dir/example/Homo_sapiens.GRCh38.cdna.21.fa

# gunzip -c ../C_long.fa.gz > C_long.fa 
# gunzip -c ../C_short.fa.gz > C_short.fa 
# gunzip -c ../GRCh38_genes_pos.bed.gz > GRCh38_genes_pos.bed 
# gunzip -c ../GRCh38_strand_pos.bed.gz > GRCh38_strand_pos.bed 
# gunzip -c ../GRCh38.21.gpd.gz > GRCh38.21.gpd 


# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test short-read alignment (HISAT2)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Index genome (chromosome 21) with HISAT2"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 hisat2-build \
# 				/work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
# 				/work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.HISAT2
# for sample in A1 A2 B1 B2
# do
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py align \
# 				--align_idx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.HISAT2 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				--1 /work_dir/${sample}_1.fq.gz  \
# 				--2 /work_dir/${sample}_2.fq.gz \
# 				--sample ${sample}
# done
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test short-read transcriptome reconstruction (StringTie)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# for sample in A1 A2 B1 B2
# do
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py reconstruct \
# 				--alignment_bam /work_dir/example/work/hisat2/${sample}/alignments.sorted.bam \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				--sample ${sample}
# done

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test quantification (Salmon-SMEM)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Index transcriptome with Salmon-SMEM"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 salmon index \
# 				-t /work_dir/example/Homo_sapiens.GRCh38.cdna.21.fa  \
# 				-i /work_dir/example/Homo_sapiens.GRCh38.cdna.21.Salmon.fmd \
# 				--type fmd
# for sample in A1 A2 B1 B2
# do
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py quantify \
# 				--quantifier_idx /work_dir/example/Homo_sapiens.GRCh38.cdna.21.Salmon.fmd \
# 				--1 /work_dir/${sample}_1.fq.gz \
# 				--2 /work_dir/${sample}_2.fq.gz \
# 				--libtype IU \
# 				--salmon_k 19 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--threads 10 \
# 				--sample ${sample} \
# 				--unzip
# done

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test differential expression analysis based on alignment"
# echo "results(DESeq2)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py diff \
# 				--alignments /work_dir/example/work/hisat2/A1/alignments.sorted.bam,/work_dir/example/work/hisat2/A2/alignments.sorted.bam /work_dir/example/work/hisat2/B1/alignments.sorted.bam,/work_dir/example/work/hisat2/B2/alignments.sorted.bam \
# 				--sample A1,A2 B1,B2 \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				--outdir /work_dir/example/out/diff-alignment/ \
# 				--workdir /work_dir/example/work/diff-alignment/

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test differential expression analysis based on alignment-free"
# echo "quantifications (DESeq2)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py diff \
# 				--quant_files /work_dir/example/work/salmon_smem/A1/quant.sf,/work_dir/example/work/salmon_smem/A2/quant.sf /work_dir/example/work/salmon_smem/B1/quant.sf,/work_dir/example/work/salmon_smem/B2/quant.sf \
# 				--sample A1,A2 B1,B2 \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				--outdir /work_dir/example/out/diff-quant \
# 				--workdir /work_dir/example/work/diff-quant

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test de novo assembly (Oases)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py denovo \
# 				--1 /work_dir/A1_1.fq.gz \
# 				--2 /work_dir/A1_2.fq.gz \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--threads 4 \
# 				--sample A1 \
# 				--file_format fastq.gz


# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test short-read fusion detection (FusionCatcher)"
# echo "Note: downloading FusionCatcher necessary data may take 1-2 hours"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# mkdir fusioncatcher_data
# cd fusioncatcher_data
# wget https://sourceforge.net/projects/fusioncatcher/files/data/human_v95.tar.gz.aa
# wget https://sourceforge.net/projects/fusioncatcher/files/data/human_v95.tar.gz.ab
# wget https://sourceforge.net/projects/fusioncatcher/files/data/human_v95.tar.gz.ac
# wget https://sourceforge.net/projects/fusioncatcher/files/data/human_v95.tar.gz.ad 
# wget https://sourceforge.net/projects/fusioncatcher/files/data/human_v95.tar.gz.ae 
# cat human_v95.tar.gz.* | tar xz
# rm -rf human_v95.tar.gz.a*

# wget http://sourceforge.net/projects/fusioncatcher/files/test/reads_1.fq.gz 
# wget http://sourceforge.net/projects/fusioncatcher/files/test/reads_2.fq.gz
# cd ..

# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py fusion \
# 				--data_dir /work_dir/example/fusioncatcher_data/human_v95/ \
# 				--input /work_dir/example/fusioncatcher_data/reads_1.fq.gz,/work_dir/example/fusioncatcher_data/reads_2.fq.gz \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--threads 4 \
# 				--sample D

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test long-read error correction (LoRDEC)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"

# docker run -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py long_correct \
# 				--kmer 23 \
# 				--solid 3 \
# 				--short /work_dir/example/C_short.fa \
# 				--long /work_dir/example/C_long.fa \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--sample C

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test long-read alignment (STARlong)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Index genome with STAR"
# mkdir STAR_genome_index_21/
# chmod -R 777 STAR_genome_index_21/
# docker run -v=${PWD}/../:/work_dir/ rnacocktail:0.3 STAR \
# 				--runMode genomeGenerate \
# 				--genomeDir /work_dir/example/STAR_genome_index_21/ \
# 				--genomeFastaFiles /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
# 				--runThreadN 4
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py long_align \
# 				--long /work_dir/example/work/lordec/C/long_corrected.fa \
# 				--threads 4 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--sample C \
# 				--genome_dir /work_dir/example/STAR_genome_index_21/

# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test long-read transcriptome reconstruction (IDP)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Index genome with STAR"
# echo "Download refFlat.txt annotation file"
# wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz
# gunzip refFlat.txt.gz
# echo "Restrict refFlat.txt annotation file to chromosome 21"
# less refFlat.txt |grep chr21 |sed 's/chr21/21/g' > refFlat.21.txt
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py align \
# 				--align_idx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.HISAT2 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf \
# 				--1 /work_dir/C_short_1.fq.gz  \
# 				--2 /work_dir/C_short_2.fq.gz \
# 				--sample C
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py long_reconstruct \
# 				--alignment /work_dir/example/work/hisat2/C/alignments.sorted.bam \
# 				--short_junction /work_dir/example/work/hisat2/C/splicesites.bed \
# 				--long_alignment /work_dir/example/work/starlong/C/Aligned.out.psl \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--threads 4 \
# 				--sample C \
# 				--ref_genome /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
# 				--ref_all_gpd /work_dir/example/GRCh38.21.gpd \
# 				--ref_gpd /work_dir/example/refFlat.21.txt \
# 				--read_length 101 \
# 				--idp /opt/IDP/src/main/python/runIDP.py


# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test long-read fusion detection (IDP-fusion)"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# wget https://www.healthcare.uiowa.edu/labs/au/IDP-fusion/files/MCF7_chr17q23-25chr20q13_Example.zip
# unzip MCF7_chr17q23-25chr20q13_Example.zip
# rm MCF7_chr17q23-25chr20q13_Example.zip
# chmod 777 -R MCF7_chr17q23-25chr20q13_Example
# echo "Download reference annotation GTF file for GRCh37"
# wget ftp://ftp.ensembl.org/pub/release-75//gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
# echo "Unzip reference annotation GTF file"
# gunzip Homo_sapiens.GRCh37.75.gtf.gz
# echo "Restrict GTF to chromosome 17 and 20"
# less Homo_sapiens.GRCh37.75.gtf |awk '{if ($1==20 || $1==17) print}'|sed 's/^/chr/' > Homo_sapiens.GRCh37.75.chromosome.17_20.gtf
# echo "Index genome (chromosome 17 and 20) with HISAT2"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 hisat2-build \
# 				/work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/genome.chr17chr20.fasta \
# 				/work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/genome.chr17chr20.HISAT2
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py align \
# 				--align_idx /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/genome.chr17chr20.HISAT2 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh37.75.chromosome.17_20.gtf \
# 				--U /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/short_reads.chr17q23-25chr20q13.tenth.fasta \
# 				--sample F \
# 				--hisat2_opts \"-f\"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py long_fusion \
# 				--alignment /work_dir/example/work/hisat2/F/alignments.sorted.bam \
# 				--short_junction /work_dir/example/work/hisat2/F/splicesites.bed \
# 				--long_alignment /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/MCF7_ready_for_fusion_20140928.chr17q23-25chr20q13.sorted.psl \
# 				--short_fasta /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/short_reads.chr17q23-25chr20q13.tenth.fasta \
# 				--long_fasta /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/long_reads.chr17q23-25chr20q13.fasta \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--threads 14 \
# 				--sample F \
# 				--ref_genome /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/genome.chr17chr20.fasta \
# 				--ref_all_gpd /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/hg19.chr17chr20.gene_est.refFlat.txt \
# 				--ref_gpd /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/refFlat_20140611.chr17chr20.gpd.txt \
# 				--read_length 89 \
# 				--genome_bowtie2_idx /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/bowtie2/bowtie2_index.chr17chr20 \
# 				--transcriptome_bowtie2_idx /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/bowtie2/txn.chr17chr20 \
# 				--uniqueness_bedgraph /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/wgEncodeDukeMapabilityUniqueness35bp.chr17chr20.bedGraph \
# 				--gmap_idx /work_dir/example/MCF7_chr17q23-25chr20q13_Example/Data/gmap_index.chr17chr20/ \
# 				--idpfusion /opt/IDP-fusion_1.1.1/bin/runIDP.py


# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Test variant calling (GATK)"
# echo "Note: the test on variant calling and RNA editing detection may take several hours"
# echo "--------------------------------------------------------"
# echo "--------------------------------------------------------"
# echo "Download reference genome FASTA file"
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# echo "Download dbsnp files"
# wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz
# gunzip All_20180418.vcf.gz
 # awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' All_20180418.vcf |sed  "s/chrMT/chrM/g" > All_20180418_chr.vcf
 # mv All_20180418_chr.vcf All_20180418.vcf
# echo "Index reference genome FASTA file"
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 samtools faidx /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.fa
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 java -jar /usr/local/bin/picard.jar CreateSequenceDictionary \
# 				R= /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# 				O= /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.dict
# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 hisat2-build \
# 				/work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# 				/work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.HISAT2
# echo "Download NA12878 FASTQ files"

# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR896/SRR896663/SRR896663_1.fastq.gz 
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR896/SRR896663/SRR896663_2.fastq.gz 

# docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py align \
# 				--align_idx /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.HISAT2 \
# 				--outdir /work_dir/example/out \
# 				--workdir /work_dir/example/work \
# 				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.gtf \
# 				--1 /work_dir/example/SRR896663_1.fastq.gz  \
# 				--2 /work_dir/example/SRR896663_2.fastq.gz \
# 				--sample E \
# 				--threads 10
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py variant \
				--alignment /work_dir/example/work/hisat2/E/alignments.sorted.bam \
				--outdir /work_dir/example/out \
				--workdir /work_dir/example/work \
				--picard /usr/local/bin/picard.jar \
				--gatk /opt/gatk-4.1.4.0/gatk-package-4.1.4.0-local.jar \
				--threads 10 \
				--sample E \
				--ref_genome /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.fa \
				--CleanSam \
				--knownsites /work_dir/example/All_20180418.vcf \
				--picard /usr/local/bin/picard.jar

echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
echo "Test RNA editing detection (GIREMI)"
echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py editing \
				--alignment /work_dir/example/work/gatk/E/bsqr.bam \
				--variant /work_dir/example/work/gatk/E/variants_filtered.vcf \
				--strand_pos /work_dir/example/GRCh38_strand_pos.bed \
				--genes_pos /work_dir/example/GRCh38_genes_pos.bed \
				--outdir /work_dir/example/out \
				--workdir /work_dir/example/work \
				--giremi_dir /usr/local/bin/ \
				--gatk /opt/gatk-4.1.4.0/gatk-package-4.1.4.0-spark.jar \
				--htslib_dir /opt/htslib-1.9/ \
				--threads 1 \
				--sample E \
				--ref_genome /work_dir/example/GRCh38_full_analysis_set_plus_decoy_hla.fa \
				--knownsites /work_dir/example/All_20180418.vcf

echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
echo "Test run all pipeline (Short-read example)"
echo "This example should successfully run for short-read alignment," 
echo "reconstruction, quantification, differential expression,"
echo "denovo assembly, variant calling, and fusion detection."
echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
cat <(less All_20180418.vcf |head -10000|grep "#") <(less All_20180418.vcf |awk '{if ($1=="chr21") print}') |sed "s/chr//g" > variants_21.vcf
cat GRCh38_genes_pos.bed |sed "s/chrM/MT/g"|sed "s/chr//g" > GRCh38_genes_pos_.bed
cat GRCh38_strand_pos.bed |sed "s/chrM/MT/g"|sed "s/chr//g" > GRCh38_strand_pos_.bed
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 samtools faidx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 java -jar /usr/local/bin/picard.jar CreateSequenceDictionary \
				R= /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
				O= /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.dict
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py all \
				--outdir /work_dir/example/out \
				--workdir /work_dir/example/work \
				--threads 10 \
				--1 /work_dir/A1_1.fq.gz,/work_dir/A2_1.fq.gz /work_dir/B1_1.fq.gz,/work_dir/B2_1.fq.gz \
				--2 /work_dir/A1_2.fq.gz,/work_dir/A2_2.fq.gz /work_dir/B1_2.fq.gz,/work_dir/B2_2.fq.gz \
				--sample all_A1,all_A2 all_B1,all_B2  \
				--ref_gtf /work_dir/example/Homo_sapiens.GRCh38.90.chromosome.21.gtf  \
				--ref_genome /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
				--align_idx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.HISAT2  \
				--quantifier_idx /work_dir/example/Homo_sapiens.GRCh38.cdna.21.Salmon.fmd \
				--unzip \
				--file_format fastq.gz \
				--CleanSam \
				--knownsites /work_dir/example/variants_21.vcf \
				--strand_pos /work_dir/example/GRCh38_strand_pos_.bed \
				--genes_pos /work_dir/example/GRCh38_genes_pos_.bed \
				--data_dir /work_dir/example/fusioncatcher_data/human_v95/ \
				--giremi_dir /usr/local/bin/ \
				--gatk /opt/gatk-4.1.4.0/gatk-package-4.1.4.0-spark.jar \
				--htslib_dir /opt/htslib-1.9/ \
				--picard /usr/local/bin/picard.jar

echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
echo "Test run all pipeline (Long-read example)"
echo "Thid example should successfully run for short-read alignment,"
echo "reconstruction, denovo assembly, quantification, "
echo "variant calling and long-read error correction," 
echo "alignment, and reconstruction"
echo "--------------------------------------------------------"
echo "--------------------------------------------------------"
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 gmap_build \
				-d /work_dir/example/gmap_chromosome.21 \
				/work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 bowtie2-build \
				/work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
				/work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 bowtie2-build \
				/work_dir/example/Homo_sapiens.GRCh38.cdna.21.fa \
				/work_dir/example/Homo_sapiens.GRCh38.cdna.21
docker run -u $UID -v=${PWD}/../:/work_dir/ rnacocktail:0.3 run_rnacocktail.py all \
				--outdir /work_dir/example/out \
				--workdir /work_dir/example/work \
				--threads 10 \
				--data_dir /work_dir/example/fusioncatcher_data/human_v95/ \
				--U /work_dir/example/C_short.fa \
				--long /work_dir/example/C_long.fa \
				--sample all_C \
				--ref_genome /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
				--star_genome_dir /work_dir/example/STAR_genome_index_21/ \
				--align_idx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21.HISAT2 \
				--quantifier_idx /work_dir/example/Homo_sapiens.GRCh38.cdna.21.Salmon.fmd \
				--gmap_idx /work_dir/example/gmap_chromosome.21 \
				--genome_bowtie2_idx /work_dir/example/Homo_sapiens.GRCh38.dna.chromosome.21 \
				--transcriptome_bowtie2_idx /work_dir/example/Homo_sapiens.GRCh38.cdna.21 \
				--ref_all_gpd /work_dir/example/GRCh38.21.gpd \
				--ref_gpd /work_dir/example/refFlat.21.txt \
				--file_format fasta \
				--CleanSam \
				--knownsites /work_dir/example/variants_21.vcf \
				--strand_pos /work_dir/example/GRCh38_strand_pos_.bed \
				--genes_pos /work_dir/example/GRCh38_genes_pos_.bed \
				--read_length 101 \
				--hisat2_opts \"-f\" \
				--giremi_dir /usr/local/bin/ \
				--gatk /opt/gatk-4.1.4.0/gatk-package-4.1.4.0-spark.jar \
				--htslib_dir /opt/htslib-1.9/ \
				--picard /usr/local/bin/picard.jar \
				--idp /opt/IDP/src/main/python/runIDP.py \
				--idpfusion /opt/IDP-fusion_1.1.1/bin/runIDP.py
