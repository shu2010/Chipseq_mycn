#!/bin/bash

# conda install -c bioconda bwa-mem2
# conda install -c bioconda trimgalore

# get genomes:
# wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
# gunzip *
# cat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna GRCh38_latest_genomic.fna > droso_human_concat_genome.fa

# conda activate

# bwa-mem2 index droso_human_concat_genome.fa 

samplelist=samplelist.txt
#readdir="/projects/drc/mycn_chip_seq_2020/"

cat $samplelist | parallel -j 8 \
	trim_galore --paired -o trim_galore/ /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed.fastq.gz \
		   /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed.fastq.gz

wait

cat $samplelist | parallel -j 48 bwa-mem2 mem -t 24 ./genomes/droso_human_concat_genome.fa \
	./trim_galore/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed_val_1.fastq.gz \
	./trim_galore/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed_val_2.fastq.gz | 
	samtools sort -o mapped/{}_nov2020.bam

wait

cat $samplelist | parallel -j 8 samtools markdup ./mapped/{}_nov2020.bam ./mapped/{}_markdup_nov2020.bam 
