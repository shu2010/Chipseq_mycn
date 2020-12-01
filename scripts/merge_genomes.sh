#!/bin/bash

# from Brundle workflow: tag and merge

# genome
sed "s/^>/>hs_/" hg38.fa > tmp
sed "s/^>/>dm_/" dm6.fa > tmp2
cat tmp tmp2 > dmhs.fa
rm tmp tmp2

# gtf
sed "s/^/hs_/" hg38.ensGene.gtf  > tmp
sed "s/^/dm_/" dm6.ensGene.gtf > tmp2
cat tmp tmp2 > dmhs.gtf
rm tmp tmp2

# blacklists
sed "s/^/hs_/" hg38-blacklist.v2.bed  > tmp
sed "s/^/dm_/" dm6-blacklist.v2.bed > tmp2
cat tmp tmp2 > dmhs-blacklist.bed
rm tmp tmp2
