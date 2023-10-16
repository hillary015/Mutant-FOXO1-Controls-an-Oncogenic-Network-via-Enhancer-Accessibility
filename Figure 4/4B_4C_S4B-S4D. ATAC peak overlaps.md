# 4B_4C_S4B-S4D. ATAC peak overlaps
The following code is an example of the code used to overlap ATAC-seq peaks.

## 4B Overlap significantly changed peaks in each cell line.
```
mergePeaks -d given 20230623_LY1_downpeaks_merge.txt 20230623_NUDUL1_downpeaks_merge.txt 20230822_DHL4_ATAC_down_overlap_1h_2h.txt > 20230822_LY1_NUD_DHL4_downpeaks_merge.txt -venn 20230822_LY1_NUD_DHL4_downpeaks_merge_venn.txt
                      ##Max distance to merge: direct overlap required (-d given)
                      ##Merging peaks...
                      ##Comparing 20230623_LY1_downpeaks_merge.txt (2302 total) and 20230623_LY1_downpeaks_merge.txt (2302 total)
                      ##Comparing 20230623_LY1_downpeaks_merge.txt (2302 total) and 20230623_NUDUL1_downpeaks_merge.txt (1162 total)
                      ##Comparing 20230623_LY1_downpeaks_merge.txt (2302 total) and 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total)
                      ##Comparing 20230623_NUDUL1_downpeaks_merge.txt (1162 total) and 20230623_LY1_downpeaks_merge.txt (2302 total)
                      ##Comparing 20230623_NUDUL1_downpeaks_merge.txt (1162 total) and 20230623_NUDUL1_downpeaks_merge.txt (1162 total)
                      ##Comparing 20230623_NUDUL1_downpeaks_merge.txt (1162 total) and 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total)
                      ##Comparing 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total) and 20230623_LY1_downpeaks_merge.txt (2302 total)
                      ##Comparing 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total) and 20230623_NUDUL1_downpeaks_merge.txt (1162 total)
                      ##Comparing 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total) and 20230822_DHL4_ATAC_down_overlap_1h_2h.txt (66 total)

mergePeaks -d given 20230623_LY1_uppeaks_merge.txt 20230623_NUDUL1_uppeaks_merge.txt 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt > 20230822_LY1_NUD_DHL4_uppeaks_merge.txt -venn 20230822_LY1_NUD_DHL4_uppeaks_merge_venn.txt
                              ##Max distance to merge: direct overlap required (-d given)
                              ##Merging peaks...
                              ##Comparing 20230623_LY1_uppeaks_merge.txt (329 total) and 20230623_LY1_uppeaks_merge.txt (329 total)
                              ##Comparing 20230623_LY1_uppeaks_merge.txt (329 total) and 20230623_NUDUL1_uppeaks_merge.txt (199 total)
                              ##Comparing 20230623_LY1_uppeaks_merge.txt (329 total) and 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total)
                              ##Comparing 20230623_NUDUL1_uppeaks_merge.txt (199 total) and 20230623_LY1_uppeaks_merge.txt (329 total)
                              ##Comparing 20230623_NUDUL1_uppeaks_merge.txt (199 total) and 20230623_NUDUL1_uppeaks_merge.txt (199 total)
                              ##Comparing 20230623_NUDUL1_uppeaks_merge.txt (199 total) and 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total)
                              ##Comparing 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total) and 20230623_LY1_uppeaks_merge.txt (329 total)
                              ##Comparing 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total) and 20230623_NUDUL1_uppeaks_merge.txt (199 total)
                              ##Comparing 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total) and 20230822_DHL4_ATAC_0v2_deseqnorm_Uppeaks_bed.txt (1 total)
```

## S4B Overlap all identified peaks 
```
cd /Volumes/Hillary_X6/FOXO1_ATAC_analysis/

mergePeaks -d given 20230623_LY1_ATAC_deseq_allpeaks_bed.txt 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart > 20230822_LY1_NUD_DHL4_allpeaks_merge.txt -venn 20230822_LY1_NUD_DHL4_allpeaks_merge_venn.txt
##Max distance to merge: direct overlap required (-d given)
        ##Merging peaks...
        ##Comparing 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total) and 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total)
        ##Comparing 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total) and 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total)
        ##Comparing 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total) and 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total)
        ##Comparing 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total) and 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total)
        ##Comparing 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total) and 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total)
        ##Comparing 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total) and 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total)
        ##Comparing 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total) and 20230623_LY1_ATAC_deseq_allpeaks_bed.txt (72344 total)
        ##Comparing 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total) and 20230623_NUDUL1_ATAC_deseqnorm_allpeaks_bed.txt (69057 total)
        ##Comparing 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total) and 20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart (107842 total)
```
## S4C Overlap significantly changed peaks in each time point by cell line
```
mergePeaks -d given 20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt 20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt > 20230822_DHL4_ATAC_down_overlap_1h_2h.txt -venn 20230822_DHL4_ATAC_down_overlap_1h_2h_venn.txt
                ##Max distance to merge: direct overlap required (-d given)
                ##Merging peaks...
                ##Comparing 20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt (21 total) and 20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt (21 total)
                ##Comparing 20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt (21 total) and 20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt (51 total)
                ##Comparing 20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt (51 total) and 20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt (21 total)
                ##Comparing 20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt (51 total) and 20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt (51 total)
```
