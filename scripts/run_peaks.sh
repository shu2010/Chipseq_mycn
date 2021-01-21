#!/bin/bash


# tried a variety of thresholds: permissive here (p=0.1) seems to provide good inupt to IDR
# however, analysis seems pretty robust to these choices 
macs2 callpeak -t blacklist_filtered/ATCACG_bl_filtered.bam -c blacklist_filtered/ACAGTG_bl_filtered.bam -f BAMPE -n MYCN_1 -g hs -p 0.1 &
macs2 callpeak -t blacklist_filtered/CGATGT_bl_filtered.bam -c blacklist_filtered/GCCAAT_bl_filtered.bam -f BAMPE -n MYCN_2 -g hs -p 0.1 &
macs2 callpeak -t blacklist_filtered/TTAGGC_bl_filtered.bam -c blacklist_filtered/CAGATC_bl_filtered.bam -f BAMPE -n MYCN_DOX_1 -g hs -p 0.1 &
macs2 callpeak -t blacklist_filtered/TGACCA_bl_filtered.bam -c blacklist_filtered/ACTTGA_bl_filtered.bam -f BAMPE -n MYCN_DOX_2 -g hs -p 0.1
