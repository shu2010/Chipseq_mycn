#!/bin/bash

samplelist=samplelist.txt

# trim reads
cat $samplelist | parallel -j 16 \
	trim_galore --paired -o trim_galore/ /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed.fastq.gz \
		   /projects/drc/mycn_chip_seq_2020/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed.fastq.gz

wait

# align with BWA-MEM2;  
cat $samplelist | parallel -j 2 "bwa-mem2 mem -t 24 ./genomes/dmhs.fa \
	./trim_galore/HCHTNCCX2_5_1_{}_150bp_271639.concat_chastity_passed_val_1.fq.gz \
	./trim_galore/HCHTNCCX2_5_2_{}_150bp_271639.concat_chastity_passed_val_2.fq.gz | 
	samtools sort -n -o mapped/{}_nov2020.bam"


wait

# fixmate for markdup
cat $samplelist | parallel -j 8 samtools fixmate -m  mapped/{}_nov2020.ns.bam mapped/{}_nov2020.fm.bam

wait

# coord sort
cat $samplelist | parallel -j 8 samtools sort -o mapped/{}_nov2020.ps.bam mapped/{}_nov2020.fm.bam

wait

# markdup
cat $samplelist | parallel -j 8 samtools markdup mapped/{}_nov2020.ps.bam mapped/{}_nov2020.markdup.bam

wait

# filter blacklist
cat $samplelist | parallel -j 8 "bedtools intersect -v -abam mapped/{}_nov2020.markdup.bam -b genomes/dmhs-blacklist.bed > blacklist_filtered/{}_bl_filtered.bam"

wait

# index - final version before peak calling
cat $samplelist | samtools index blacklist_filtered/{}_bl_filtered.bam
