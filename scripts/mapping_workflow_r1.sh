#!/bin/bash

# conda activate

# conda install -c bioconda bwa-mem2
# conda install -c bioconda trimgalore

# get genomes:

# mkdir -p genomes
# cd genomes

# wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
# gunzip *
# cat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna GRCh38_latest_genomic.fna > droso_human_concat_genome.fa
# bwa-mem2 index droso_human_concat_genome.fa 

# cd ..

# mkdir -p mapped
# mkdir -p markdup
# mkdir -p flagstats

# ls /projects/drc/mycn_chip_seq_2020/*.gz | sed 's/^.*5_1_\|^.*5_2_\|_150.*$//g' | sort -u > samplelist.txt

samplelist=samplelist.txt
# readdir="/projects/drc/mycn_chip_seq_2020/"

cat $samplelist | parallel -j 8 \
	trim_galore --paired -o trim_galore/ /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed.fastq.gz \
		   /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed.fastq.gz

wait

cat $samplelist | parallel -j 2 "bwa-mem2 mem -t 24 ./genomes/droso_human_concat_genome.fa \
	./trim_galore/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed_val_1.fq.gz \
	./trim_galore/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed_val_2.fq.gz | 
	samtools sort -o mapped/{}_nov2020.bam"

wait

cat $samplelist | parallel -j 8 samtools flagstat ./mapped/{}_nov2020.bam > ./flagstats/{}_nov2020.bam.flagstat
cat $samplelist | parallel -j 8 samtools markdup ./mapped/{}_nov2020.bam ./markdup/{}_markdup_nov2020.bam 
