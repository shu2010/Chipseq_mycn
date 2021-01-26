# take from MACS2 github issue #356: how to incorporate control
# specify read length
read_len=150

# specify peaks cut off
log_qval_cutoff=1.301
log_pval_cutoff=2

# remove duplicates
macs2 filterdup -f BAMPE -i ${prefix}.sorted.bam --keep-dup=1 -o ${prefix}_filterdup.bed

# convert sam to bam for control sample
samtools view -S -b $ctrl_file | samtools sort - -o ${ctrl_prefix}.sorted.bam
# remove duplicates
macs2 filterdup -f BAMPE -i ${ctrl_prefix}.sorted.bam --keep-dup=1 
    -o ${ctrl_prefix}_filterdup.bed

## generate ChIP coverage tracks
# get d (fragment size) from each peak file then calculate p normalization factor
# for background
d=$(cat $peaks_dir/${prefix}_peaks.xls | awk '/# d/ {print $4}')
d2=$(echo $d/2 | bc)
macs2 pileup -f BEDPE -B -i ${prefix}_filterdup.bed 
     -o ${prefix}_filterdup.pileup.bdg --extsize $d

## build local background track from control
# d background
macs2 pileup -f BEDPE -B -i ${ctrl_prefix}_filterdup.bed 
    -o ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg --extsize $d2
# slocal background
macs2 pileup -f BEDPE -B --extsize 500 -i ${ctrl_prefix}_filterdup.bed 
     -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg
# normalize the small local (slocal) background
p=$(echo "scale=4;"$d/1000 | bc)
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_1k_bg.pileup.bdg -m multiply -p $p 
    -o ${ctrl_prefix}_${prefix}_filterdup_1k_bg_norm.bdg
# large local background default 10k (--extsize 5000)
macs2 pileup -f BEDPE -i ${ctrl_prefix}_filterdup.bed -B --extsize 5000 
    -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg.bdg
# normalize the large local (llocal) background
p=$(echo "scale=4;"$d/10000 | bc)
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_10k_bg.bdg -m multiply -p $p 
    -o ${ctrl_prefix}_${prefix}_filterdup_10k_bg_norm.bdg
# combine and generate maximum background noise
macs2 bdgcmp -m max -t ${ctrl_prefix}_${prefix}_filterdup_1k_bg_norm.bdg 
    -c ${ctrl_prefix}_${prefix}_filterdup_10k_bg_norm.bdg 
    -o ${ctrl_prefix}_${prefix}_filterdup_1k_10k_bg_norm.bdg
# Then, take the maximum then by comparing with d background:
macs2 bdgcmp -m max -t ${ctrl_prefix}_${prefix}_filterdup_1k_10k_bg_norm.bdg 
    -c ${ctrl_prefix}_${prefix}_filterdup_d_bg.pileup.bdg 
    -o ${ctrl_prefix}_${prefix}_filterdup_d_1k_10k_bg_norm.bdg

## finally combine with genome wide background
#p=the_number_of_control_reads*fragment_length/genome_size
ctrl_num_read=$(wc -l ${ctrl_prefix}_filterdup.bed | awk '{print $1}') 
genome_size=2700000000
p=$(echo "scale=4;"$ctrl_num_read*$d/$genome_size | bc)
echo $ctrl_num_read $genome_size $p
macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_d_1k_10k_bg_norm.bdg -m max -p $p
    -o ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg		

# scaling chip and control based on spike in
treat_num_read=$(wc -l ${prefix}_filterdup.bed | awk '{print $1}')
treat_spike_num_read=$(wc -l ${prefix}.ecoli_filterdup.bed | awk '{print $1}')
ctrl_spike_num_read=$(wc -l ${ctrl_prefix}.ecoli_filterdup.bed | awk '{print $1}')
min_spike_num_read=$(wc -l *ecoli_filterdup.bed | sort -nk 1 | head -n 1 | awk '{print $1}')
scale_factor_spike_treat=$(echo "scale=4;"$min_spike_num_read/$treat_spike_num_read | bc)
scale_factor_spike_ctrl=$(echo "scale=4;"$min_spike_num_read/$ctrl_spike_num_read | bc)
p=$(echo "scale=4;"$scale_factor_spike_treat/$scale_factor_spike_ctrl | bc)
echo spike in scale factor=$p
if (( $(echo $p '>' 1 |bc -l) )); then
	p=$(echo "scale=10;"1/$p | bc)
	echo scale down control...
	# scale down control
	$run macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg -m multiply -p $p \
            -o ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg
	$run macs2 bdgopt -i ${prefix}_filterdup.pileup.bdg -m multiply -p 1 -o ${prefix}_filterdup_scale.pileup.bdg
	echo scaled spike-in control= $(echo $p*$ctrl_spike_num_read | bc)
	echo scaled spike-in treat= $(echo $treat_spike_num_read | bc)
else
	echo scale down treatment
	# scale down treatment
	$run macs2 bdgopt -i ${prefix}_filterdup.pileup.bdg -m multiply -p $p -o ${prefix}_filterdup_scale.pileup.bdg
	$run macs2 bdgopt -i ${ctrl_prefix}_${prefix}_filterdup_local_bias_raw.bdg -m multiply -p 1 \
            -o ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg
	echo scaled spike-in control= $(echo $p*$treat_spike_num_read | bc)
	echo scaled spike-in treat= $(echo $ctrl_spike_num_read | bc)
fi

# compare chip and control and calculate q and pvalue
macs2 bdgcmp -t ${prefix}_filterdup_scale.pileup.bdg 
       -c  ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg -m qpois 
        -o ${prefix}_qvalue.bdg
macs2 bdgcmp -t ${prefix}_filterdup_scale.pileup.bdg 
       -c  ${ctrl_prefix}_${prefix}_filterdup_local_lambda.bdg -m ppois 
       -o ${prefix}_pvalue.bdg
#Call peaks on score track using a cutoff
macs2 bdgpeakcall -i ${prefix}_qvalue.bdg -c $log_qval_cutoff -l $d -g $read_len 
        -o ${prefix}_qval_peaks.bed
macs2 bdgpeakcall -i ${prefix}_pvalue.bdg -c $log_pval_cutoff -l $d -g $read_len 
        -o ${prefix}_pval_peaks.bed
`

edited: to fix mixed up scaling.
p.s: with `-f BEDPE`, you don't really need to specify `--extsize` anymore.
