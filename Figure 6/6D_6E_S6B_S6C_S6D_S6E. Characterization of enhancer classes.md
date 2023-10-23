# 6D_S6B_S6C_S6D_S6E. Characterization of enhancer classes
The following code is an example of the code used to identify and classify the enhancer classes in bash and R.
## Identification of enhancer classes
Functional peaks, ATAC peaks downregulated in the DHL4s, and all ATAC-seq peaks identified in the DHL4s were intersected using bedtools in bash.
```
#Overlap functional peaks with downregulated DHL4 ATAC peaks
bedtools intersect -a 20230703_funcFKHmotifs_nohead.bed -b ../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt ../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt -u -filenames > 20230822_intersect_functionalpeaks_withboundmotifs_wsigDHL4peaks_u_overlaponly.bed
bedtools intersect -a ../20230703_funcFKHmotifs_nohead.bed -b ../../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_0v1_deseqnorm_downpeaks_bed.txt ../../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_0v2_deseqnorm_downpeaks_bed.txt -v > ../202308223_intersect_functionalpeaks_withsigDHL4peaks_v_nooverlap.bed
#Overlap functional peaks without downregulated DHL4 peak with all DHL4 ATAC peaks
bedtools intersect -a ../202308223_intersect_functionalpeaks_withsigDHL4peaks_v_nooverlap.bed -b ../../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart -v > ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_vnooverlap.txt
bedtools intersect -a ../202308223_intersect_functionalpeaks_withsigDHL4peaks_v_nooverlap.bed -b ../../../FOXO1_ATAC_analysis/20230822_DHL4_ATAC_deseqnorm_allpeaks_bed.txt.filepart -u > ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt
```
## Determine which motifs are DIV2 vs. FKH
DIV2 motifs were identified using bedtools in bash.
```
cat motifs_LY1_0.25hrvs0h/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2/beds/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2_all.bed motifs_NUDUL1_0.25hvs0h/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2/beds/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2_all.bed motifs_DHL4_0.5hvs0h/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2/beds/TRANSFAC.FOXO1_DIV2_FOXO1_DIV2_all.bed > 20230828_DIV2_motifs_all_cat.txt
sort -k1,1 -k2,2n 20230828_DIV2_motifs_all_cat.txt > 20230828_DIV2_motifs_all_cat_sorted.bed
bedtools merge -i 20230828_DIV2_motifs_all_cat_sorted_bed4.bed > 20230828_DIV2_motifs_all_cat_sorted_merged.bed

bedtools intersect -a 20230801_funcFKHmotifs_motifonly.bed -b 20230828_DIV2_motifs_all_cat_sorted_merged.bed -wb > 20230828_intersect_functionalmotifs_wDIV2all_wb.txt
bedtools intersect -a 20230703_funcFKHmotifs_nohead.bed -b 20230828_DIV2_motifs_all_cat_sorted_merged.bed -wa -wb > 20230828_intersect_functionalpeaks_wDIV2all_wa_wb.txt
```
## Create bar graph of average # of FKH motifs per enhancer class
The average # of FKH motifs per enhancer class was plotted in R
### Load libraries
```{r}
library(tidyverse)
library(hillaryscolors)
```

### Load data
```{r}
intersected_peaks <- read_delim("20230822_intersect_functionalpeaks_withboundmotifs_wa_wb.bed", col_names = c("chr_peak", "start_peak", "end_peak", "chr_motif", "start_motif", "end_motif"))

intersected_peaks_DHL4 <- read_delim("20230822_intersect_functionalpeaks_withboundmotifs_wsigDHL4peaks_wa_wb.bed", col_names = c("chr_peak", "start_peak", "end_peak", "file", "chr_motif", "start_motif", "end_motif"))

DIV2_motifs <- read_delim("20230828_intersect_functionalpeaks_wDIV2all_wa_wb.txt", col_names = c("chr_peak", "start_peak", "end_peak", "chr_motif", "start_motif", "end_motif"))
```
### Rearrange data
```{r}
intersected_peaks_DHL4_adj <- intersected_peaks_DHL4 %>%
  mutate(functional_peak = paste(chr_peak, start_peak, end_peak, sep = "_")) %>%
  filter(file != "20230801_funcFKHmotifs_motifonly.bed")

DIV2_motifs_adj <- DIV2_motifs %>%
    mutate(functional_motif = paste(chr_motif, start_motif, end_motif, sep = "_"))

intersected_peaks_adj <- intersected_peaks %>%
  mutate(functional_peak = paste(chr_peak, start_peak, end_peak, sep = "_")) %>%
  mutate(functional_motif = paste(chr_motif, start_motif, end_motif, sep = "_")) %>%
  mutate(motifs = case_when(functional_motif %in% DIV2_motifs_adj$functional_motif ~ "2",  !(functional_motif %in% DIV2_motifs_adj$functional_motif) ~ "1")) %>%
  group_by(functional_peak) %>%
  summarise(motifs_per_peak = sum(as.numeric(motifs))) %>% 
  mutate(DHL4_sig = case_when(functional_peak %in% intersected_peaks_DHL4_adj$functional_peak ~ TRUE, TRUE ~ FALSE))

pwt_motifs_per_peak <- pairwise.wilcox.test(intersected_peaks_adj$motifs_per_peak, intersected_peaks_adj$DHL4_sig,
                 p.adjust.method = "BH")
print(pwt_motifs_per_peak)

DHL4_sig <- intersected_peaks_adj %>% 
  filter(DHL4_sig == TRUE) %>% 
  select(motifs_per_peak) %>%
  mutate(motifs_per_peak = as.numeric(motifs_per_peak))

print(summary(DHL4_sig))

DHL4_nonsig <- intersected_peaks_adj %>% 
  filter(DHL4_sig == FALSE) %>% 
  select(motifs_per_peak) %>%
  mutate(motifs_per_peak = as.numeric(motifs_per_peak))

print(summary(DHL4_nonsig))

wt_motifs_per_peak <- wilcox.test(DHL4_sig$motifs_per_peak, DHL4_nonsig$motifs_per_peak, alternative = "two.sided")
print(wt_motifs_per_peak)
```
### Bar graph
```{r}
bargraph <- intersected_peaks_adj %>%
  group_by(DHL4_sig) %>%
  mutate(mean = mean(motifs_per_peak), sd = sd(motifs_per_peak), sem = sd/sqrt(n()))
print(bargraph)

bg <- ggplot(bargraph, aes(x = DHL4_sig, y = mean)) + 
  geom_bar(stat ="identity", position = position_dodge(0.65), color = "black", width = 0.5, fill = "#AFAFAF") +
  geom_point(aes(x = DHL4_sig, y = motifs_per_peak, fill = DHL4_sig), stat ="identity", position = position_jitterdodge(jitter.height = 0.2, dodge.width = 0.65), color = "black", ) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sem), color = "Black", width = .2, position = position_dodge(0.65)) +
  ggtitle("Bound Motifs Per Peak") +
  labs(x = "Peak Set", y = "# of Motifs") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 15, by = 5), limits = c(0, 15)) +
  scale_x_discrete(limits = pos) +
  theme_hillary() +
  theme(plot.title = element_text(vjust = 2)) +
  coord_cartesian(clip = "off")
bg

dotplot  <- ggplot(bargraph, aes(x = DHL4_sig, y = motifs_per_peak, group = DHL4_sig)) + 
  geom_point(size = 1, position = position_jitter(width = 0.1)) +
  ggtitle("Motifs per peak") +
  labs(x = "Class", y = "Motifs per peak") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 15),  breaks = seq(0, 15, 3)) +
  scale_x_discrete(limits = pos) +
  geom_errorbar(aes(ymax = mean, ymin = mean), colour = "black", width = .5) +
    theme_hillary() +
  scale_x_discrete(limits = pos) +
  theme(plot.title = element_text(vjust = 2)) +
  coord_cartesian(clip = "off")
dotplot
```

### Print to png
```{r}
Cairo::Cairo(file = "20230822_boundmotifsperfunctionalpeak_hist_4x8.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 8, 
             height = 4, 
             pointsize = 12, 
             dpi = 300)
hist
dev.off()

Cairo::Cairo(file = "20230822_boundmotifsperfunctionalpeak_boxplot_4x4.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 4, 
             height = 4, 
             pointsize = 12, 
             dpi = 300)
boxplot
dev.off()
```
## BINDetect (WT vs. Mutant)
Differentially bound motifs between mut. and WT cell lines were identified using TOBIAS in bash.
```
TOBIAS BINDetect --motifs  ${PEAK_DIR}/coremotifs_jaspar_plusDIV2.txt --signals /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/DHL4_ATAC_0h_footprints.bw /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/LY1_FOXO1_0hr_merged.sorted_corrected_footprints.bw.filepart --genome ${GENOME_DIR}/hg19/hg19.fa --peaks ${PEAK_DIR}/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed --outdir ${PROCESSED_DIR}/Functional_peaks/DHL40hvLY10h --cond_names DHL4 LY1 --cores 8 --verbosity 1
TOBIAS BINDetect --motifs  ${PEAK_DIR}/coremotifs_jaspar_plusDIV2.txt --signals /Volumes/Hillary_X6/Nud_DHL4_ATAC_forTobias/TF_footprinting_out/DHL4_ATAC_0h_footprints.bw ${PEAK_DIR}/TF_footprinting_out/NUDUL1_ATAC_0h_footprints.bw --genome ${GENOME_DIR}/hg19/hg19.fa --peaks ${PEAK_DIR}/TF_footprinting_out/20230703_funcFKHmotifs_nohead.bed --outdir ${PROCESSED_DIR}/Functional_peaks/DHL40hvNUD0h --cond_names DHL4 NUD --cores 8 --verbosity 1
```
## Volcano plot
Volcano plots of differentially bound motifs were plotted in R.
## Load Libraries
```{r}
library(tidyverse)
library(hillaryscolors)
library(stringi)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
```

## Load bindetect results from TOBIAS and expressed genes
```{r}
bindetect_DHL4_vNUD <- readr::read_delim("bindetect_results_DHL4vNUD_0h.txt")
bindetect_DHL4_vLY1 <- readr::read_delim("bindetect_results_DHL4_vLY1_0h.txt")

expressed_genes_NUDUL1 <- readr::read_delim("../../Exp 0460 NUDUL1_FOXO1_PS/20230205 NUDUL1_FOXO1_PS_expressedgenes_duprm_manfilter.txt")
expressed_genes_LY1 <- readr::read_delim("../../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/20220912 LY1_FOXO1_PS_expressedgenes_duprm_manfilter.txt")
expressed_genes_DHL4 <- readr::read_delim("../../Exp 0071/20230530 DHL4_FOXO1_PS_expressedgenes_duprm_manfilter.txt")
```

### Select genes expressed in all three cell lines
```{r}
expressed_all <- expressed_genes_LY1 %>%
  full_join(expressed_genes_NUDUL1, by = "Transcript") %>%
  full_join(expressed_genes_DHL4, by = "Transcript") %>%
  mutate(Gene_all = coalesce(Gene.x, Gene.y, Gene))
```

### Adjust file 
```{r}
tf_results_adj <- function(df, changecol, pcol) {
  df %>%
  dplyr::mutate(change = changecol * -1) %>%
  dplyr::mutate(pval = pcol) %>%
  dplyr::mutate(motif = str_extract(name, "((?<=\\.[:digit:].)[[:print:]]+)")) %>%
  dplyr::mutate(motif = replace_na(motif, "DIV2")) %>%
  dplyr::filter(motif %in% expressed_all$Gene_all | motif == "DIV2") %>%
  dplyr::mutate(label = case_when(str_detect(motif, "FOX") ~ "FKH", str_detect(motif, "DIV2") ~ "DIV2", str_detect(motif, "POU") ~ "POU", str_detect(motif, "MEF") ~ "MEF", str_detect(motif, "FOX", negate = TRUE) & str_detect(motif, "POU", negate = TRUE) & str_detect(motif, "MEF", negate = TRUE) ~ "zNA")) %>%
  dplyr::filter(motif_id != "MA0032.1" & motif_id != "MA0033.1" & motif_id != "MA0479.1")   %>%
  dplyr::mutate(top2.5percent_pval = quantile(pval, probs = 0.025)) %>%
  dplyr::mutate(top2.5percent_l2fc = quantile(change, probs = 0.025)) %>%
  dplyr::mutate(top5percent_pval = quantile(pval, probs = 0.05)) %>%
  dplyr::mutate(top5percent_l2fc = quantile(change, probs = 0.05)) %>%
  dplyr::mutate(bottom2.5percent_l2fc = quantile(change, probs = 0.975)) %>%
  dplyr::mutate(bottom5percent_l2fc = quantile(change, probs = 0.95)) %>%
  dplyr::mutate(sig_2.5 = case_when(pval <= top2.5percent_pval & change > 0 ~ "Up", pval <= top2.5percent_pval &  change < 0 ~ "Down", pval < 1 & change <= top2.5percent_l2fc ~ "Down", pval < 1 &  change >= bottom2.5percent_l2fc ~ "Up", TRUE ~ "zNS")) %>%
  dplyr::mutate(sig_5 = case_when(pval <= top5percent_pval & change > 0 ~ "Up", pval <= top5percent_pval &  change < 0 ~ "Down", pval < 1 & change <= top5percent_l2fc ~ "Down", pval < 1 &  change >= bottom5percent_l2fc ~ "Up", TRUE ~ "zNS"))
}

DHL4vNud <- tf_results_adj(bindetect_DHL4_vNUD, bindetect_DHL4_vNUD$`DHL4_NUD_change`,  bindetect_DHL4_vNUD$`DHL4_NUD_pvalue`)
DHL4vLY1 <- tf_results_adj(bindetect_DHL4_vLY1, bindetect_DHL4_vLY1$`DHL4_LY1_change`,  bindetect_DHL4_vLY1$`DHL4_LY1_pvalue`)
```
### Make volcano plot (POU FKH label color)
```{r}
volplot_selectlabel <- function(df, title, ylimhi) {
  ggplot(df %>% dplyr::arrange(desc(label))) +
  geom_point(aes(x = change, y = -log10(pval), fill = label, alpha = label, colour = label), shape = 21,  size = 4) +  
  ggtitle(title) +
  xlab("Change in Binding Score") +
  ylab("-log10(p-value)") +
  scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_hillary()+
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("#EF7A85","#9EE493","#B892FF", "#0FA3B1", "#AFAFAF")) +
  scale_color_manual(values = c("#000000", "#000000", "#000000", "#000000", "#AFAFAF")) +
  theme(plot.title = element_text(hjust = 0.5, size = "20", face = "bold")) +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, .8))
}
## Make average plot
WTvmut <- select(DHL4vLY1, motif_id, label, change, pval)%>%
  left_join(select(DHL4vNud, motif, motif_id, label, change, pval), by = "motif_id") %>%
  mutate(average_change = (change.x + change.y)/2, average_p = (pval.x + pval.y)/2)

WTvmut_volplot <- ggplot(WTvmut %>% dplyr::arrange(desc(label.x))) +
  geom_point(aes(x = average_change, y = -log10(average_p), fill = label.x, alpha = label.x, colour = label.x), shape = 21,  size = 4) +  
  ggtitle("WT vs Mut TFBS Change") +
  xlab("Change in Binding Score") +
  ylab("-log10(p-value)") +
  scale_y_continuous(limits = c(0, 125), expand = c(0,0)) +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_hillary()+
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("#EF7A85","#9EE493","#B892FF", "#0FA3B1", "#AFAFAF")) +
  scale_color_manual(values = c("#000000", "#000000", "#000000", "#000000", "#AFAFAF")) +
  theme(plot.title = element_text(hjust = 0.5, size = "20", face = "bold")) +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, .8))
WTvmut_volplot
```

### Write to PNG
```{r}
png_6x5 <- function(title, plot) {
  Cairo::Cairo(file=title, 
             bg="white",
             type="png",
             units="in", 
             width=6, 
             height=5, 
             pointsize=12, 
             dpi=300)
plot
}

png_6x5("20230814_DHL4vNud_0h_TF_footprinting_functional_peaks_volplot.png", DHL4_vNUD_vol)
dev.off()
png_6x5("20230814_DHL4vLY1_0h_TF_footprinting_functional_peaks_volplot.png", DHL4_vLY1_vol)
dev.off()
png_6x5("20230814_DHL4vNud_0h_TF_footprinting_functional_peaks_volplot_colorpalPOU_FKH_DIV2_MEF.png", DHL4_vNUD_vol_selectlabel)
dev.off()
png_6x5("20230814_DHL4vLY1_0h_TF_footprinting_functional_peaks_volplot_colorpalPOU_FKH_DIV2_MEF.png", DHL4_vLY1_vol_selectlabel)
dev.off()
png_6x5("20230814_DHL4vavgNudLY1_0h_TF_footprinting_functional_peaks_volplot_colorpalPOU_FKH_DIV2_MEF.png", WTvmut_volplot)
dev.off()
```
## Intersect class 2 enhancers with POU/MEF2 motifs
Class 2 enhancer peaks were intersected with POU/MEF2 motifs using bedtools and HOMER in bash.
```
cat DHL40hvLY10h/MA0507.1.POU2F2_MA0507.1/beds/MA0507.1.POU2F2_MA0507.1_DHL4_bound.bed DHL40hvNUD0h/MA0507.1.POU2F2_MA0507.1/beds/MA0507.1.POU2F2_MA0507.1_DHL4_bound.bed DHL40hvLY10h/MA0507.2.POU2F2_MA0507.2/beds/MA0507.2.POU2F2_MA0507.2_DHL4_bound.bed DHL40hvNUD0h/MA0507.2.POU2F2_MA0507.2/beds/MA0507.2.POU2F2_MA0507.2_DHL4_bound.bed DHL40hvLY10h/MA0628.1.POU6F1_MA0628.1/beds/MA0628.1.POU6F1_MA0628.1_DHL4_bound.bed DHL40hvNUD0h/MA0628.1.POU6F1_MA0628.1/beds/MA0628.1.POU6F1_MA0628.1_DHL4_bound.bed DHL40hvLY10h/MA0785.1.POU2F1_MA0785.1/beds/MA0785.1.POU2F1_MA0785.1_DHL4_bound.bed DHL40hvNUD0h/MA0785.1.POU2F1_MA0785.1/beds/MA0785.1.POU2F1_MA0785.1_DHL4_bound.bed  DHL40hvLY10h/MA1115.1.POU5F1_MA1115.1/beds/MA1115.1.POU5F1_MA1115.1_DHL4_bound.bed DHL40hvNUD0h/MA1115.1.POU5F1_MA1115.1/beds/MA1115.1.POU5F1_MA1115.1_DHL4_bound.bed DHL40hvLY10h/MA1549.1.POU6F1_MA1549.1/beds/MA1549.1.POU6F1_MA1549.1_DHL4_bound.bed DHL40hvNUD0h/MA1549.1.POU6F1_MA1549.1/beds/MA1549.1.POU6F1_MA1549.1_DHL4_bound.bed > 20230823_POU_DHL4_bound_motifs_cat.bed

sort -k1,1 -k2,2n 20230823_POU_DHL4_bound_motifs_cat.bed > 20230823_POU_DHL4_bound_motifs_cat_sorted.bed
bedtools merge -i 20230823_POU_DHL4_bound_motifs_cat_sorted.bed > 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed

bedtools intersect -a ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt -b 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed -wa -wb > 20230823_intersect_functionalpeaksmutonly_DHL4conpeaks_withPOUboudmoitfs_wa_wb.txt

bedtools intersect -a ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt -b 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed -wa -wb > 20230823_intersect_functionalpeaksmutonly_DHL4conpeaks_withPOUboudmoitfs_wa_wb.txt
bedtools intersect -a ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt -b 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed -u > 20230823_intersect_functionalpeaksmutonly_DHL4conpeaks_withPOUboudmoitfs_ufuncpeakonly.txt
bedtools intersect -a ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt -b 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed -wb > 20230823_intersect_functionalpeaksmutonly_DHL4conpeaks_withPOUboudmoitfs_wbmotifonly.txt

cat DHL40hvLY10h/MA0052.1.MEF2A_MA0052.1/beds/MA0052.1.MEF2A_MA0052.1_DHL4_bound.bed DHL40hvNUD0h/MA0052.1.MEF2A_MA0052.1/beds/MA0052.1.MEF2A_MA0052.1_DHL4_bound.bed DHL40hvLY10h/MA0052.2.MEF2A_MA0052.2/beds/MA0052.2.MEF2A_MA0052.2_DHL4_bound.bed DHL40hvNUD0h/MA0052.2.MEF2A_MA0052.2/beds/MA0052.2.MEF2A_MA0052.2_DHL4_bound.bed DHL40hvLY10h/MA0052.3.MEF2A_MA0052.3/beds/MA0052.3.MEF2A_MA0052.3_DHL4_bound.bed DHL40hvNUD0h/MA0052.3.MEF2A_MA0052.3/beds/MA0052.3.MEF2A_MA0052.3_DHL4_bound.bed DHL40hvLY10h/MA0052.4.MEF2A_MA0052.4/beds/MA0052.4.MEF2A_MA0052.4_DHL4_bound.bed DHL40hvNUD0h/MA0052.4.MEF2A_MA0052.4/beds/MA0052.4.MEF2A_MA0052.4_DHL4_bound.bed DHL40hvLY10h/MA0497.1.MEF2C_MA0497.1/beds/MA0497.1.MEF2C_MA0497.1_DHL4_bound.bed DHL40hvNUD0h/MA0497.1.MEF2C_MA0497.1/beds/MA0497.1.MEF2C_MA0497.1_DHL4_bound.bed DHL40hvLY10h/MA0660.1.MEF2B_MA0660.1/beds/MA0660.1.MEF2B_MA0660.1_DHL4_bound.bed DHL40hvNUD0h/MA0660.1.MEF2B_MA0660.1/beds/MA0660.1.MEF2B_MA0660.1_DHL4_bound.bed DHL40hvLY10h/MA0773.1.MEF2D_MA0773.1/beds/MA0773.1.MEF2D_MA0773.1_DHL4_bound.bed DHL40hvNUD0h/MA0773.1.MEF2D_MA0773.1/beds/MA0773.1.MEF2D_MA0773.1_DHL4_bound.bed > 20230823_MEF_DHL4_bound_motifs_cat.bed

sort -k1,1 -k2,2n 20230823_MEF_DHL4_bound_motifs_cat.bed > 20230823_MEF_DHL4_bound_motifs_cat_sorted.bed
bedtools merge -i 20230823_MEF_DHL4_bound_motifs_cat_sorted.bed > 20230823_MEF_DHL4_bound_motifs_cat_sorted_merged.bed

mergePeaks -d given ../20230823_intersect_functionalpeaksmutonly_withDHL4conpeaks_uoverlap.txt 20230823_POU_DHL4_bound_motifs_cat_sorted_merged.bed 20230823_MEF_DHL4_bound_motifs_cat_sorted_merged.bed > 20230823_mergePeaks_funcmutonlywithDHL4conpeak_POU2DHL4bound_MEF2DHL4bound.txt -venn 20230823_mergePeaks_funcmutonlywithDHL4conpeak_POU2DHL4bound_MEF2DHL4bound_venn.txt
```
## Heatmaps of POU/MEF2 factor expression
Heatmaps of POU/MEF2 factor expression were created in R from RNA-seq data.
### Load Library
```{r}
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(viridis)
library(hillaryscolors)
library(ggbreak)
```

### Load Data
```{r}
DHL4_FOXO1_0hrv6hr_geneexp_diff <- read.delim("../Exp 0131_FOXO1_ETO_RS/DHL4_FOXO1_0hrv6hr_geneexp_diff.txt")
DHL4_FOXO1_0hrv24hr_geneexp_diff <- read.delim("../Exp 0131_FOXO1_ETO_RS/DHL4_FOXO1_0hrv24hr_geneexp_diff.txt")
DHL4_FOXO1_0hrv6hr_counts <- read.delim("../Exp 0131_FOXO1_ETO_RS/DHL4_FOXO1_0hrv6hr_genes_read_group_tracking.txt") %>%
  mutate(Cell_line = "DHL4") %>%
  mutate(condition = case_when(condition == "0hr" ~ "0hr", condition == "4hr" ~ "6hr"))
DHL4_FOXO1_0hrv24hr_counts <- read.delim("../Exp 0131_FOXO1_ETO_RS/DHL4_FOXO1_0hrv24hr_genes_read_group_tracking.txt") %>%
  mutate(Cell_line = "DHL4") %>%
  mutate(condition = case_when(condition == "0hr" ~ "0hr", condition == "4hr" ~ "24hr"))

LY1_FOXO1_0hrv6hr_geneexp_diff <- read.delim("../Exp 0131_FOXO1_ETO_RS/LY1_FOXO1_0hrv6hr_geneexp_diff.txt")
LY1_FOXO1_0hrv24hr_geneexp_diff <- read.delim("../Exp 0131_FOXO1_ETO_RS/LY1_FOXO1_0hrv24hr_geneexp_diff.txt")
LY1_FOXO1_0hrv6hr_counts <- read.delim("../Exp 0131_FOXO1_ETO_RS/LY1_FOXO1_0hrv6hr_genes_read_group_tracking.txt") %>%
  mutate(Cell_line = "LY1")
LY1_FOXO1_0hrv24hr_counts <- read.delim("../Exp 0131_FOXO1_ETO_RS/LY1_FOXO1_0hrv24hr_genes_read_group_tracking.txt") %>%
  mutate(Cell_line = "LY1")

NUD_FOXO1_FKBP_0hrv6hr_geneexp_diff <- read.delim("NUD_FOXO1_FKBP_0v6h_gene_exp_diff.txt")
NUD_FOXO1_FKBP_0hrv24hr_geneexp_diff <- read.delim("NUD_FOXO1_FKBP_0v24h_gene_exp_diff.txt")
NUD_FOXO1_FKBP_0hrv6hr_counts <- read.delim("NUD_FOXO1_FKBP_0v6_genes.read_group_tracking.txt") %>%
  mutate(Cell_line = "NUD")
NUD_FOXO1_FKBP_0hrv24hr_counts <- read.delim("NUD_FOXO1_FKBP_0v24_genes.read_group_tracking.txt") %>%
  mutate(Cell_line = "NUD")
```

### Join Tables by gene_id
```{r}
POU <- LY1_FOXO1_0hrv6hr_geneexp_diff %>% 
  left_join(select(LY1_FOXO1_0hrv24hr_geneexp_diff, gene_id, LY1_T24 = value_2)) %>%
  left_join(select(DHL4_FOXO1_0hrv6hr_geneexp_diff, gene_id, DHL4_T0 = value_1, DHL4_T6 = value_2)) %>%
  left_join(select(DHL4_FOXO1_0hrv24hr_geneexp_diff, gene_id, DHL4_T24 = value_2)) %>%
  left_join(select(NUD_FOXO1_FKBP_0hrv6hr_geneexp_diff, gene_id, NUD_T0 = value_1, NUD_T6 = value_2)) %>%
  left_join(select(NUD_FOXO1_FKBP_0hrv24hr_geneexp_diff, gene_id, NUD_T24 = value_2)) %>%
  select(gene_id, NUD_T0, NUD_T6, NUD_T24, LY1_T0 = value_1, LY1_T6 = value_2, LY1_T24, DHL4_T0, DHL4_T6, DHL4_T24) %>%
  filter(str_detect(gene_id, "POU")) %>%
  filter(gene_id != "POU5FIP" & gene_id != "POU5F1P3" & gene_id != "POU5F1P4")

MEF2 <- LY1_FOXO1_0hrv6hr_geneexp_diff %>% 
  left_join(select(LY1_FOXO1_0hrv24hr_geneexp_diff, gene_id, LY1_T24 = value_2)) %>%
  left_join(select(DHL4_FOXO1_0hrv6hr_geneexp_diff, gene_id, DHL4_T0 = value_1, DHL4_T6 = value_2)) %>%
  left_join(select(DHL4_FOXO1_0hrv24hr_geneexp_diff, gene_id, DHL4_T24 = value_2)) %>%
  left_join(select(NUD_FOXO1_FKBP_0hrv6hr_geneexp_diff, gene_id, NUD_T0 = value_1, NUD_T6 = value_2)) %>%
  left_join(select(NUD_FOXO1_FKBP_0hrv24hr_geneexp_diff, gene_id, NUD_T24 = value_2)) %>%
  select(gene_id, NUD_T0, NUD_T6, NUD_T24, LY1_T0 = value_1, LY1_T6 = value_2, LY1_T24, DHL4_T0, DHL4_T6, DHL4_T24) %>%
  filter(str_detect(gene_id, "MEF2")) %>%
  filter(gene_id != "PRAMEF2" & gene_id != "PRAMEF20,PRAMEF21" & gene_id != "PRAMEF21" & gene_id != "PRAMEF22,PRAMEF3" & gene_id != "MEF2BNB" & gene_id != "MEF2BNB-MEF2B")
```

### Make heatmap with z-score or FPKM
```{r}
rnames_POU = POU [ , 1]
mat_data_POU = data.matrix(POU [ , c(2,5,8)])
rownames(mat_data_POU) = rnames_POU
my_palette = colorRamp2(c(0, 2.34375, 4.6875, 9.375, 18.75, 37.5, 75, 150, 300), c("#FFFFFF", "#EDF8B1FF", "#C7E9B4FF", "#7FCDBBFF", "#41B6C4FF", "#1D91C0FF", "#225EA8FF", "#253494FF", "#081D58FF"))

heatmap_POU <- Heatmap(mat_data_POU, col = my_palette, column_title = "POU Factors", cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 18), rect_gp = gpar(col= "black"))
heatmap_POU

rnames_MEF2 = MEF2 [ , 1]
mat_data_MEF2 = data.matrix(MEF2 [ , c(2,5,8)])
rownames(mat_data_MEF2) = rnames_MEF2
my_palette = colorRamp2(c(0, 2.34375, 4.6875, 9.375, 18.75, 37.5, 75, 150, 300), c("#FFFFFF", "#EDF8B1FF", "#C7E9B4FF", "#7FCDBBFF", "#41B6C4FF", "#1D91C0FF", "#225EA8FF", "#253494FF", "#081D58FF"))

heatmap_MEF2 <- Heatmap(mat_data_MEF2, col = my_palette, column_title = "MEF2 Factors", cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 18), rect_gp = gpar(col= "black"))
heatmap_MEF2
```

### Save as PNG
```{r}
Cairo::Cairo(file="20230817 NUD_LY1_DHL4_FOXO1_RS_POU_heatmap_5x10.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=10, 
             pointsize=12, 
             dpi=300)
heatmap_POU
dev.off()

Cairo::Cairo(file="20230823 NUD_LY1_DHL4_FOXO1_RS_MEF2_heatmap_5x10.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=10, 
             pointsize=12, 
             dpi=300)
heatmap_MEF2
dev.off()
```
