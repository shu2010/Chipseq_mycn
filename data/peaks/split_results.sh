# The BAMs and BEDs should be split into Drosophila (for normalization) and into Human (for diffbind)



### BAMs, split in two



echo `pwd`

mkdir human
mkdir drosophila


for xls in *peaks.xls
do
grep -v "dm_chr" $xls | sed s/hs_chr/chr/g > human/$xls
grep -v "hs_chr" $xls | sed s/dm_chr/chr/g > drosophila/$xls
done

for narrow in *narrowPeak
do
grep -v "dm_chr" $narrow | sed s/hs_chr/chr/g > human/$narrow
grep -v "hs_chr" $narrow | sed s/dm_chr/chr/g > drosophila/$narrow
done
