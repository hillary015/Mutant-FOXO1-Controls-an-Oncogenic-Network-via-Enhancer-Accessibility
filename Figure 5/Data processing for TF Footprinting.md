# Transcription Factor Foot Printing Data Processing
The code used to process ATAC-seq data for transcription factor footprinting with TOBIAS is presented below.

## Merge BAM files
BAM files from biological replicates were merged, sorted, and indexed using SAMtools in bash.
```
samtools merge -f 9330-NUDUL1_ATAC_0h_merged.bam 9330-HL-0001_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam 9330-HL-0002_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam
samtools sort -o DHL4_ATAC_0h_merged_sorted.bam -T temp_ -@ 10 DHL4_ATAC_0h_merged.bam
samtools index DHL4_ATAC_0h_merged_sorted.bam
```
## ATACorrect
Merged BAM files were corrected for Tn5 sequence bias using TOBIAS in bash.
```
PEAK_DIR=/Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/
PROCESSED_DIR=/Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/
GENOME_DIR=/Volumes/Hillary_X6/genomes/
#this first command is correcting for the bias of Tn5 cut sites, your output is a directory containing a *corrected.bw file and some others
if test -d ${PROCESSED_DIR}; then  echo "${PROCESSED_DIR} exists"; else mkdir ${PROCESSED_DIR} && echo "${PROCESSED_DIR} created"; fi

for i in 0h 0.25h 0.5h 1h 2h
do
TOBIAS ATACorrect --bam ${PEAK_DIR}/DHL4_ATAC_${i}_merged_sorted.bam --genome ${GENOME_DIR}/hg19/hg19.fa --peaks ${PEAK_DIR}/DHL4_consensuspeaks.bed.txt \
    --blacklist ${GENOME_DIR}/hg19_blacklist/hg19.blacklistpeaks.bed --outdir ${PROCESSED_DIR}/DHL4_ATAC_${i} --cores 8
```
## Score Bigwig
Corrected bigwig files were scored using TOBIAS in bash.
```
TOBIAS ScoreBigwig --signal ${PROCESSED_DIR}/DHL4_ATAC_${i}/DHL4_ATAC_${i}_merged_sorted_corrected.bw --regions ${PEAK_DIR}/DHL4_consensuspeaks.bed.txt \
    --output ${PROCESSED_DIR}/DHL4_ATAC_${i}_footprints.bw --cores 8
```
