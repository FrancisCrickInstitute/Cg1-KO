#!/bin/bash
# Consensus peaks (>=2 WT, >=3 HET replicates; 50% reciprocal overlap) + featureCounts.


ROOT=$(pwd)/results_atacseq/bwa/merged_library   
PEAKDIR=${ROOT}/macs2/broad_peak
OUT=${PEAKDIR}/overlaps-1
mkdir -p ${OUT}/counts_table

WT="WT_Skin_REP1 WT_Skin_REP2 WT_Skin_REP3"
HET="HET_Skin_REP1 HET_Skin_REP2 HET_Skin_REP3 HET_Skin_REP4"

bp() { for s in "$@"; do echo ${PEAKDIR}/${s}.mLb.clN_peaks.broadPeak; done; }

# --- WT: peaks present in >=2 replicates ---
cat $(bp $WT) | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT}/merged_wt.bed
bedtools intersect -a ${OUT}/merged_wt.bed -b $(bp $WT) \
  -f 0.50 -r -c -nonamecheck | awk '$NF >= 2' > ${OUT}/overlap_wt.bed

# --- HET: peaks present in >=3 replicates ---
cat $(bp $HET) | sort -k1,1 -k2,2n | bedtools merge -i - > ${OUT}/merged_het.bed
bedtools intersect -a ${OUT}/merged_het.bed -b $(bp $HET) \
  -f 0.50 -r -c -nonamecheck | awk '$NF >= 3' > ${OUT}/overlap_het.bed

# --- merge WT and HET consensus peaks ---
bedops -m ${OUT}/overlap_wt.bed ${OUT}/overlap_het.bed > ${OUT}/merged_overlap_peaks.bed

# SAF for featureCounts, BED with peak IDs for annotation
awk -F '\t' 'BEGIN{OFS=FS}{$2=$2+1; id="merged_overlap_macs2_"++n; print id,$1,$2,$3,"."}' \
  ${OUT}/merged_overlap_peaks.bed > ${OUT}/merged_overlap-1_peaks.saf
awk -F '\t' 'BEGIN{OFS=FS}{$2=$2+1; id="merged_overlap_macs2_"++n; print $1,$2,$3,id,"0","."}' \
  ${OUT}/merged_overlap_peaks.bed > ${OUT}/merged_overlap-1_peaksid.bed

# --- featureCounts over consensus peaks (order: HET then WT) ---
BAMS=""
for s in $HET $WT; do BAMS="$BAMS ${ROOT}/${s}.mLb.clN.sorted.bam"; done

featureCounts -F SAF -O --fracOverlap 0.2 -p -T 16 -s 0 \
  -a ${OUT}/merged_overlap-1_peaks.saf \
  -o ${OUT}/counts_table/merged_peaks_macs2_overlap-1.counts $BAMS

# drop the leading '#' comment line so R reads it cleanly
awk 'NR>1' ${OUT}/counts_table/merged_peaks_macs2_overlap-1.counts \
  > ${OUT}/counts_table/merged_peaks_macs2_overlap-1.counts.tsv
