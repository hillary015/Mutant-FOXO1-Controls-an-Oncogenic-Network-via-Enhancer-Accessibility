#S4F. Histone modification overlaps
The following code is an example of the code used to overlap histone modification ChIP-seq/CUT&RUN peaks with ATAC-seq peaks.

## S4F Overlap histone modification peaks from control samples with significantly downregulated ATAC peaks.
Peaks were overlapped using HOMER in bash.
```
mergePeaks -d given 20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt  /Volumes/Hillary_X6/5437-LY1_H3K4me3_0.1xlinked/20240117\ LY1_FOXO1_H3K4me3_consensuspeaks_0v6_mo1_bed.txt /Volumes/Hillary_X6/5400-LY1_H3K27ac/20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_consensuspeaks_mo1_bed.txt /Volumes/Hillary_2p/Hatzi_NatImmuno_2019_LY1_BCL6_LSD1_H3K4me1_ChIP/LY1_H3K4me1_NT.peaks.bed > 20240119_mergePeaks_downATAC_H3K4me3_mo1_H3K27ac_mo1_H3K4me1.txt -venn 20240119_mergePeaks_downATAC_H3K4me3_mo1_H3K27ac_mo1_H3K4me1_venn.txt
```
