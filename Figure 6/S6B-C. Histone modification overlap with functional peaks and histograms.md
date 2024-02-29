#S6B-C. Histone modification overlap with functional peaks and histograms
The following code is an example of the code used to the histone modifications with functional peaks and create histogrmas of histone modification signal around functional peaks.

## 6B Overlap significantly changed peaks in each cell line.
Histone modification peaks from the control sample were overlapped with functional peaks using HOMER in bash.
```
mergePeaks -d given /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed /Volumes/Hillary_X6/5437-LY1_H3K4me3_0.1xlinked/20240117\ LY1_FOXO1_H3K4me3_consensuspks_0v6_mo1_bed.txt /Volumes/Hillary_X6/5400-LY1_H3K27ac/20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_consensuspeaks_mo1_bed.txt /Volumes/Hillary_2p/Hatzi_NatImmuno_2019_LY1_BCL6_LSD1_H3K4me1_ChIP/LY1_H3K4me1_NT.peaks.bed > 20240114_mergePeaks_funcFKH_H3K4me3_mo1_H3K27ac_mo1_H3K4me1.txt -venn 20240117_mergePeaks_funcFKH_H3K4me3_mo1_H3K27ac_mo1_H3K4me1_venn.txt
```

## S4C Histograms
### Identify motifs with specific histone modifications
Functional motifs were intersected with individual histone modifications using Bedtools in bash.
```
bedtools intersect -a /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed -b /Volumes/Hillary_X6/5437-LY1_H3K4me3_0.1xlinked/20240109_LY1_FOXO1_H3K4me3_peaks_bed.txt -u -wa -wb >  20240109_funcFKH_overlap_H3K4me3_overlap_u.txt
bedtools intersect -a /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed -b /Volumes/Hillary_2p/Hatzi_NatImmuno_2019_LY1_BCL6_LSD1_H3K4me1_ChIP/LY1_H3K4me1_NT.peaks.bed -u -wa -wb >  20240109_funcFKH_overlap_H3K4me1_overlap_u.txt
bedtools intersect -a /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed -b /Volumes/Hillary_X6/LY1_FOXO1_H3K27ac_0_2hr_comb_exps/20240104_0164_0264_LY1_FOXO1_H3K27ac_peaks.bed -u -wa -wb >  20240109_funcFKH_overlap_H3K27ac_overlap_u.txt
```

## Generate histogram data
Histogram data was generated using deepTools in conda.
```
## Functional Peaks
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230801_funcFKHmotifs_motifonly.bed -S 5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw  --skipZeros --binSize 10 -o 20240114_H3K3K4me3_0v6_computematrix_funcFKHpeaks_ba3000_bs10.gz

plotProfile --matrixFile 20240114_H3K3K4me3_0v6_computematrix_funcFKHpeaks_ba3000_bs10.gz --outFileName 20240114_H3K4me3_0v6_functionalFKH_bs10.svg --plotHeight 15 --plotWidth 15 --plotType lines --perGroup --yAxisLabel "Signal" --regionsLabel "Functional FKH Motifs" --refPointLabel "Motif Center" --outFileNameData 20240114_H3K4me3_0v6_computematrix_funcFKHpeaks_ba3000_bs10_profile.txt

## Functional Peaks with histone modification
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R /Volumes/Hillary_X6/histonemark_analysis/20240117_funcFKH_motifs_overlap_H3K4me3_overlap_u.txt -S 5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw  --skipZeros --binSize 10 -o 20240117_H3K4me3_0v6_computematrix_funcFKHpeaks_H3K4me3only_ba3000_bs10.gz

plotProfile --matrixFile 20240114_H3K4me3_0v6_computematrix_funcFKHmotifs_H3K4me3only_ba3000_bs10.gz --outFileName 20240114_H3K4me3_0v6_functionalFKHmotiffs_H3K4me3_bs10.svg --plotHeight 15 --plotWidth 15 --plotType lines --perGroup --yAxisLabel "Signal" --regionsLabel "Functional FKH Motifs" --refPointLabel "Motif Center" --outFileNameData 20240117_H3K4me3_0v6_computematrix_funcFKHmotifs_H3K4me3only_ba3000_bs10_profile.txt

## Non-functional Peaks
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230801_nonfuncFKHmotifs_motifonly.bed -S 5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw 5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10.bw  --skipZeros --binSize 10 -o 20240117_H3K4me3_0v6_computematrix_nonfuncFKHmotifs_ba3000_bs10.gz

plotProfile --matrixFile 20240117_H3K4me3_0v6_computematrix_nonfuncFKHmotifs_ba3000_bs10.gz --outFileName 20240117_H3K4me3_0v6_computematrix_nonfuncFKHmotifs_ba3000_bs10_profile.svg --plotHeight 15 --plotWidth 15 --plotType lines --perGroup --yAxisLabel "Signal" --regionsLabel "Functional FKH Motifs" --refPointLabel "Motif Center" --outFileNameData 20240117_H3K4me3_0v6_computematrix_nonfuncFKHmotifs_ba3000_bs10_profile.txt
```

### Make histograms in R
Files were adjusted in Excel so that column names were populated correctly in R.
```
## Load libraries
library(tidyverse)
library(hillaryscolors)
library(viridis)

## Load Files
funcFKH <- read.delim("20240114_H3K4me3_0v6_computematrix_funcFKHpeaks_ba3000_bs10_profile_adj.txt") %>%
  select(1:602) %>%
  pivot_longer(X1:X600, names_to = "bin", values_to = "score") %>%
  mutate(bin = parse_number(bin)) %>%
  pivot_wider(names_from = Samples, values_from = score) %>%
  mutate(avg0 = (`5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2) %>%
  mutate(avg6 = (`5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2)

funcFKH_H3K4me3only <- read.delim("20240117_H3K4me3_0v6_computematrix_funcFKHmotifs_H3K4me3only_ba3000_bs10_profile_adj.txt") %>%
  select(1:602) %>%
  pivot_longer(X1:X600, names_to = "bin", values_to = "score") %>%
  mutate(bin = parse_number(bin)) %>%
  pivot_wider(names_from = Sample, values_from = score) %>%
  mutate(avg0 = (`5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2) %>%
  mutate(avg6 = (`5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2)

nonfuncFKH <- read.delim("20240117_H3K4me3_0v6_computematrix_nonfuncFKHmotifs_ba3000_bs10_profile_adj.txt") %>%
  select(1:602) %>%
  pivot_longer(X1:X600, names_to = "bin", values_to = "score") %>%
  mutate(bin = parse_number(bin)) %>%
  pivot_wider(names_from = Samples, values_from = score) %>%
  mutate(avg0 = (`5437-HL-1_LY1_FOXO1_H3K4me3_0hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-2_LY1_FOXO1_H3K4me3_0hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2) %>%
  mutate(avg6 = (`5437-HL-5_LY1_FOXO1_H3K4me3_6hrA_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10` + `5437-HL-6_LY1_FOXO1_H3K4me3_6hrB_xlinked.hg19sc.F4q10.sorted.scrm.sfdeseqrecip.bs10`)/2)

## Make histograms
histpalette <- c("0 hr" = "#325D58", "6 hr" = "#539792")

funcFKH_hist <- ggplot (funcFKH, aes(x = bin), color = ) +
  geom_line(aes(y = avg0, color = "0 hr"), size = 0.75) +
  geom_line(aes(y = avg6, color = "6 hr"), size = 0.75) +
  scale_x_continuous(breaks = seq(1, 601, by = 300), expand = c(0, 0), limits = c(1, 601)) +
  scale_y_continuous(breaks = seq(0, 20, by = 10), expand = c(0, 0), limits = c(0, 20)) +
  ggtitle("H3K4me3 all Func") +
  labs(x = "Distance from Motif Center", y = "H3K4me3 Signal") +
  theme_hillary() +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold")) +
  scale_color_manual(values = histpalette)
funcFKH_hist

funcFKH_H3K4me3_hist <- ggplot (funcFKH_H3K4me3only, aes(x = bin), color = ) +
  geom_line(aes(y = avg0, color = "0 hr"), size = 0.75) +
  geom_line(aes(y = avg6, color = "6 hr"), size = 0.75) +
  scale_x_continuous(breaks = seq(1, 601, by = 300), expand = c(0, 0), limits = c(1, 601)) +
  scale_y_continuous(breaks = seq(0, 20, by = 10), expand = c(0, 0), limits = c(0, 20)) +
  ggtitle("H3K4me3 Func w/H3K4me3") +
  labs(x = "Distance from Motif Center", y = "H3K4me3 Signal") +
  theme_hillary() +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold")) +
  scale_color_manual(values = histpalette)
funcFKH_H3K4me3_hist

nonfuncFKH_hist <- ggplot (nonfuncFKH, aes(x = bin), color = ) +
  geom_line(aes(y = avg0, color = "0 hr"), size = 0.75) +
  geom_line(aes(y = avg6, color = "6 hr"), size = 0.75) +
  scale_x_continuous(breaks = seq(1, 601, by = 300), expand = c(0, 0), limits = c(1, 601)) +
  scale_y_continuous(breaks = seq(0, 20, by = 10), expand = c(0, 0), limits = c(0, 20)) +
  ggtitle("H3K4me3 nonFunc") +
  labs(x = "Distance from Motif Center", y = "H3K4me3 Signal") +
  theme_hillary() +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold")) +
  scale_color_manual(values = histpalette)
nonfuncFKH_hist

## Save to png
png <- function(name, data) {
ggplot2::ggsave(file= name, 
             device="png",
             scale = 1,
             dpi=300,
             width = 5,
             height = 4,
             plot = data)
}

smallpng <- function(name, data) {
ggplot2::ggsave(file= name, 
             device="png",
             scale = 1,
             dpi=300,
             width = 3.5,
             height = 2.5,
             plot = data)
}

png("20240117_histo_H3K4me3sig_funcpeaks_bin5_3kb_5x4.png", funcFKH_hist)
png("20240117_histo_H3K4me3sig_funcpeakswH3K4me3_bin5_3kb_5x4.png", funcFKH_H3K4me3_hist)
png("20240117_histo_H3K4me3sig_nonfuncpeaks_bin5_3kb_5x4.png", nonfuncFKH_hist)

smallpng("20240117_histo_H3K4me3sig_funcpeaks_bin5_3kb_3.5x2.5.png", funcFKH_hist)
smallpng("20240117_histo_H3K4me3sig_funcpeakswH3K4me3_bin5_3kb_3.5x2.5.png", funcFKH_H3K4me3_hist)
smallpng("20240117_histo_H3K4me3sig_nonfuncpeaks_bin5_3kb_3.5x2.5.png", nonfuncFKH_hist)
```
