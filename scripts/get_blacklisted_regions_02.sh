#!/bin/bash
##get blacklisted regions to exclude from BAMS

wget -c https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip ENCFF356LFX.bed.gz
mv ENCFF356LFX.bed hg38_blacklist.bed

# Drosophila dm3
#wget -c https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz
#gunzip dm6-blacklist.v2.bed.gz

### Generate a merged Drosophila/Human blacklist
droso=dm6-blacklist.v2.bed
homo=hg38_blacklist.bed

sed "s/^/hs_/" $homo |cut -f1-3 > tmp
sed "s/^/dm_/" $droso > tmp2
cat tmp tmp2 > dmhs-blacklist.bed
rm tmp tmp2
