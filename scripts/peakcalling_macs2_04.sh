### MACS2 peak caller

### Run macs on the blacklisted data
##please assign appropriate bams for treatment and control
##using the sample sheet will be useful
##Example from Alex data
##for PE data
#conda activate macs2

#macs2 callpeak -t HCHTNCCX2_5_TTAGGC.sorted.marked_duplicates.bam –-outdir peak_dir -c HCHTNCCX2_5_CAGATC.sorted.marked_duplicates.bam -f BAMPE -g hs --bdg -n dox_Rep1


mkdir ./peaks
cd peaks
control=../blacklist_filtered/input_DNA.bam

for bam in ../blacklist_filtered/*bam
do
root=`basename $bam .bam`

macs2 callpeak -t $bam -c $control -f BAMPE -n $root -g hs --bdg
done
