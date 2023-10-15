# 2D_S2D
The following code was used to generate pausing index heatmaps and bargraphs from PRO-seq data in R.
## Load Libraries
```{r}
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(xlsx)
library(viridis)
library(hillaryscolors)
```

## Load files
```{r}
LY1_PI_0v0.5_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_0v0.5.txt") 
LY1_PI_0v1_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_0v1.txt") 
LY1_PI_2_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_0v2.txt") 

DHL4_PI_0v0.5_Deseq <- readr::read_delim("../Exp 0071/DHL4_pindex_0v0.5.txt") 
DHL4_PI_0v1_Deseq <- readr::read_delim("../Exp 0071/DHL4_pindex_0v1.txt") 
DHL4_PI_2_Deseq <- readr::read_delim("../Exp 0071/pindex_0v2.txt") 

NUD_PI_0v1_Deseq <- readr::read_delim("pindex_0v1h.txt") 
NUD_PI_0v2_Deseq <- readr::read_delim("pindex_0v2.txt")

LY1_PI_change_0v0.5_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_change_0v0.5.txt") 
LY1_PI_change_0v1_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_change_0v1.txt") 
LY1_PI_change_2_Deseq <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/pindex_change_0v2.txt") 

DHL4_PI_change_0v0.5_Deseq <- readr::read_delim("../Exp 0071/DHL4_pindex_change_0v0.5.txt") 
DHL4_PI_change_0v1_Deseq <- readr::read_delim("../Exp 0071/DHL4_pindex_change_0v1.txt") 
DHL4_PI_change_2_Deseq <- readr::read_delim("../Exp 0071/DHL4_pindex_change_0v2.txt") 

NUD_PI_change_0v1_Deseq <- readr::read_delim("pindex_change_0v1h.txt") 
NUD_PI_change_0v2_Deseq <- readr::read_delim("pindex_change_0v2.txt")

clusters <- readr::read_delim("20230606 PS_gbs_heatmap_NUD_LY1_DHL4_siggenes_LY1_NUD_overlap_replicatesgbd_km3_clusters.txt") %>%
  dplyr::select(Transcript, Gene, cluster) %>%
  dplyr::distinct()
```

## Combine pi tbles
```{r}
km3clustering_PI <- clusters %>%
  dplyr::left_join(LY1_PI_0v0.5_Deseq, by = "Transcript") %>%
  dplyr::left_join(LY1_PI_0v1_Deseq, by = "Transcript") %>%
  dplyr::left_join(LY1_PI_2_Deseq, by = "Transcript") %>%
  dplyr::left_join(NUD_PI_0v1_Deseq, by = "Transcript") %>%
  dplyr::left_join(NUD_PI_0v2_Deseq, by = "Transcript") %>%
  dplyr::left_join(NUD_PI_0v2_Deseq, by = "Transcript") %>%
  dplyr::left_join(DHL4_PI_0v0.5_Deseq, by = "Transcript") %>%
  dplyr::left_join(DHL4_PI_0v1_Deseq, by = "Transcript") %>%
  dplyr::left_join(DHL4_PI_2_Deseq, by = "Transcript") %>%
  dplyr::select(cluster, Transcript, NUD0A = `9322-HL-0001_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.y`, NUD0B =  `9322-HL-0002_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.y`, NUD1A = `9322-HL-0004_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, NUD1B = `9322-HL-0005_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, NUD2A = `9322-HL-0006_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.y`, NUD2B = `9322-HL-0007_S1_L001_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.y`, LY10A = `7771-HL-1_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.x`, LY10B =  `7771-HL-2_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex.x`, LY10.5A = `7771-HL-3_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, LY10.5B =  `7771-HL-4_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, LY11A = `7771-HL-5_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, LY11B = `7771-HL-6_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, LY12A = `7771-HL-7_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, LY12B = `7771-HL-8_comb.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted-pindex`, DHL40A = `3466-KS-2-DHL4_FOXO1_FKBP_0hr-1.hg19spikeinrm2.F4q10.sorted-pindex.x`, DHL40B = `3466-KS-3-DHL4_FOXO1_FKBP_0hr-2.hg19spikeinrm2.F4q10.sorted-pindex.x`, DHL40.5A = `3466-KS-8-DHL4_FOXO1_FKBP_30min-1.hg19spikeinrm2.F4q10.sorted-pindex`, DHL40.5B = `3466-KS-9-DHL4_FOXO1_FKBP_30min-2.hg19spikeinrm2.F4q10.sorted-pindex`, DHL41A = `3466-KS-4-DHL4_FOXO1_FKBP_1hr-1.hg19spikeinrm2.F4q10.sorted-pindex`, DHL41B = `3466-KS-5-DHL4_FOXO1_FKBP_1hr-2.hg19spikeinrm2.F4q10.sorted-pindex`, DHL42A = `3466-KS-10-DHL4_FOXO1_FKBP_2hr-1.hg19spikeinrm2.F4q10.sorted-pindex`, DHL42B = `3466-KS-11-DHL4_FOXO1_FKBP_2hr-2.hg19spikeinrm2.F4q10.sorted-pindex`) %>%
  dplyr::mutate(l2fcNUD1A = log2((NUD1A/((NUD0A + NUD0B)/2)))) %>%
  dplyr::mutate(l2fcNUD1B = log2((NUD1B/((NUD0A + NUD0B)/2)))) %>%
  dplyr::mutate(l2fcNUD2A = log2((NUD2A/((NUD0A + NUD0B)/2)))) %>%
  dplyr::mutate(l2fcNUD2B = log2((NUD2B/((NUD0A + NUD0B)/2)))) %>%
  dplyr::mutate(l2fcLY10.5A = log2((LY10.5A/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcLY10.5B = log2((LY10.5B/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcLY11A = log2((LY11A/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcLY11B = log2((LY11B/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcLY12A = log2((LY12A/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcLY12B = log2((LY12B/((LY10A + LY10B)/2)))) %>%
  dplyr::mutate(l2fcDHL40.5A = log2((DHL40.5A/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate(l2fcDHL40.5B = log2((DHL40.5B/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate(l2fcDHL41A = log2((DHL41A/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate(l2fcDHL41B = log2((DHL41B/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate(l2fcDHL42A = log2((DHL42A/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate(l2fcDHL42B = log2((DHL42B/((DHL40A + DHL40B)/2)))) %>%
  dplyr::mutate_all(~replace(., is.infinite(.), 0))
```

## Make Heatmaps
```{r, fig.show = 'hold', fig.dim = c(2, 4)}
#All heatmaps
my_palette_pi = colorRamp2(c(-2, 0, 2), c("#053061", "#f5f5f5", "#70000C"))

pindex_siggenes <-as.data.frame(km3clustering_PI)

rnames_pi = pindex_siggenes [ , 2]
mat_data_pi = data.frame(pindex_siggenes [ , 25:40])
rownames(mat_data_pi) = rnames_pi
cluster = pindex_siggenes [ , 1]

heatmap_pi <- Heatmap(mat_data_pi, col = my_palette_pi, column_title = "LY1 FOXO1 Pausing Index, siggenes any tp km3 clusters", cluster_columns = FALSE, cluster_rows = FALSE, row_split = cluster, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 4),  row_gap = unit(2, "mm"))
heatmap_pi
```

## Print to PNG
```{r, fig.show = 'hide'}
Cairo::Cairo(file="202300906_NUD_LY1_DHL4_km3_FOXO1_PROseqsiggenes_pausingindexchange_l2fc_10x5.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=10, 
             pointsize=12, 
             dpi=300)
heatmap_pi
dev.off()
```

## Combine pi change tables
```{r}
pindex_change <- dplyr::left_join(clusters, LY1_PI_change_0v0.5_Deseq, by = "Transcript", suffix = c("_km3", "_0.5")) %>%
 dplyr::left_join(LY1_PI_change_0v1_Deseq, by = "Transcript", suffix = c("_0.5", "_1")) %>%
 dplyr::left_join(LY1_PI_change_2_Deseq, by = "Transcript", suffix = c("_1", "_2")) %>%
 dplyr::left_join(NUD_PI_change_0v1_Deseq, by = "Transcript", suffix = c("_2", "_N1")) %>%
 dplyr::left_join(NUD_PI_change_0v2_Deseq, by = "Transcript", suffix = c("_N1", "_N2")) %>%
 dplyr::left_join(DHL4_PI_change_0v0.5_Deseq, by = "Transcript", suffix = c("_N2", "_D0.5")) %>%
 dplyr::left_join(DHL4_PI_change_0v1_Deseq, by = "Transcript", suffix = c("_D0.5", "_D1")) %>%
 dplyr::left_join(DHL4_PI_change_2_Deseq, by = "Transcript", suffix = c("_D1", "_D2")) %>%
 dplyr::select(Gene_km3, cluster, Transcript, log2fc_0.5, log2fc_1, log2fc_2, log2fc_N1, log2fc_N2, log2fc_N1, log2fc_D0.5, log2fc_D1, log2fc_D2, FDR_0.5, FDR_1, FDR_2, FDR_N1, FDR_N2, FDR_D0.5, FDR_D1, FDR_D2) %>%
 dplyr::mutate_all(~replace(., is.na(.), 0))
```

## ID sig changes in pausing index per cluster
```{r}
siggenes_wpiclass <- pindex_change %>%
   dplyr::mutate(pi_0.5hrsig = case_when(FDR_0.5 >= 0.05 ~ "NS", log2fc_0.5 > -0.585 & log2fc_0.5 < 0.585 ~ "NS", FDR_0.5 < 0.05 & log2fc_0.5 > 0.585 ~ "Up", FDR_0.5 < 0.05 & log2fc_0.5 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_1hrsig = case_when(FDR_1 >= 0.05 ~ "NS", log2fc_1 > -0.585 & log2fc_1 < 0.585 ~ "NS", FDR_1 < 0.05 & log2fc_1 > 0.585 ~ "Up", FDR_1 < 0.05 & log2fc_1 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_2hrsig = case_when(FDR_2 >= 0.05 ~ "NS", log2fc_2 > -0.585 & log2fc_2 < 0.585 ~ "NS", FDR_2 < 0.05 & log2fc_2 > 0.585 ~ "Up", FDR_2 < 0.05 & log2fc_2 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_N1hrsig = case_when(FDR_N1 >= 0.05 ~ "NS", log2fc_N1 > -0.585 & log2fc_N1 < 0.585 ~ "NS", FDR_N1 < 0.05 & log2fc_N1 > 0.585 ~ "Up", FDR_N1 < 0.05 & log2fc_N1 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_N2hrsig = case_when(FDR_N2 >= 0.05 ~ "NS", log2fc_N2 > -0.585 & log2fc_N2 < 0.585 ~ "NS", FDR_N2 < 0.05 & log2fc_N2 > 0.585 ~ "Up", FDR_N2 < 0.05 & log2fc_N2 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_D0.5hrsig = case_when(FDR_D0.5 >= 0.05 ~ "NS", log2fc_D0.5 > -0.585 & log2fc_D0.5 < 0.585 ~ "NS", FDR_D0.5 < 0.05 & log2fc_D0.5 > 0.585 ~ "Up", FDR_D0.5 < 0.05 & log2fc_D0.5 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_D1hrsig = case_when(FDR_D1 >= 0.05 ~ "NS", log2fc_D1 > -0.585 & log2fc_D1 < 0.585 ~ "NS", FDR_D1 < 0.05 & log2fc_D1 > 0.585 ~ "Up", FDR_D1 < 0.05 & log2fc_D1 < -0.585 ~ "Down")) %>%
   dplyr::mutate(pi_D2hrsig = case_when(FDR_D2 >= 0.05 ~ "NS", log2fc_D2 > -0.585 & log2fc_D2 < 0.585 ~ "NS", FDR_D2 < 0.05 & log2fc_D2 > 0.585 ~ "Up", FDR_D2 < 0.05 & log2fc_D2 < -0.585 ~ "Down"))

siggenes_wpichange_cluster_count_0.5 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_0.5hrsig = replace_na(pi_0.5hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_0.5hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_0.5hrsig, values_from = genes) %>%
  dplyr::mutate(Down = 0) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_1 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_1hrsig = replace_na(pi_1hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_1hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_1hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_2 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_2hrsig = replace_na(pi_2hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_2hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_2hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = 0) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_N1 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_N1hrsig = replace_na(pi_N1hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_N1hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_N1hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_N2 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_N2hrsig = replace_na(pi_N2hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_N2hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_N2hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_D0.5 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_D0.5hrsig = replace_na(pi_D0.5hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_D0.5hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_D0.5hrsig, values_from = genes) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_D1 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_D1hrsig = replace_na(pi_D1hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_D1hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_D1hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_cluster_count_D2 <- siggenes_wpiclass %>%
  dplyr::mutate(pi_D2hrsig = replace_na(pi_D2hrsig, "NS")) %>%
  dplyr::group_by(cluster, pi_D2hrsig) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = pi_D2hrsig, values_from = genes) %>%
  dplyr::mutate(Up = replace_na(Up, 0)) %>%
  dplyr::mutate(Down = replace_na(Down, 0)) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_nocluster <- siggenes_wpiclass %>%
  mutate(NUD = case_when(pi_N1hrsig == "Down" | pi_N2hrsig == "Down" ~ "Down", pi_N1hrsig == "Up" | pi_N2hrsig == "Up" ~ "Up",  pi_N1hrsig == "NS" | pi_N2hrsig == "NS" ~ "NS")) %>%
  mutate(LY1 = case_when(pi_0.5hrsig == "Down" | pi_1hrsig == "Down" | pi_2hrsig == "Down" ~ "Down", pi_0.5hrsig == "Up" | pi_1hrsig == "Up" | pi_2hrsig == "Up" ~ "Up", pi_0.5hrsig == "NS"  | pi_1hrsig == "NS" | pi_2hrsig == "NS" ~ "NS")) %>%
  mutate(DHL4 = case_when(pi_D0.5hrsig == "Down" | pi_D1hrsig == "Down" | pi_D2hrsig == "Down" ~ "Down", pi_D0.5hrsig == "Up" | pi_D1hrsig == "Up" | pi_D2hrsig == "Up" ~ "Up", pi_D0.5hrsig == "NS"  | pi_D1hrsig == "NS" | pi_D2hrsig == "NS" ~ "NS"))  %>%
pivot_longer(NUD:DHL4, names_to = "cell_line", values_to = "group") %>%
  dplyr::group_by(cell_line, group) %>%
  dplyr::summarise(genes = n()) %>%
  tidyr::pivot_wider(names_from = group, values_from = genes) %>%
  dplyr::mutate(totalgenespercluster = (Down + Up + NS)) %>%
  dplyr::mutate(C_percentdown = (Down/totalgenespercluster)*100) %>%
  dplyr::mutate(B_percentup = (Up/totalgenespercluster)*100) %>%
  dplyr::mutate(A_percentns = (NS/totalgenespercluster)*100) %>%
  tidyr::pivot_longer(cols = C_percentdown:A_percentns, names_to = "group", values_to = "percent")

siggenes_wpichange_gbc <- siggenes_wpiclass %>%
  mutate(gbc_class = case_when(cluster == 1 | cluster == 2 ~ "gbd", cluster == 3 ~ "gbu")) %>%
  mutate(NUD = case_when(pi_N1hrsig == "Down" | pi_N2hrsig == "Down" ~ "Down", pi_N1hrsig == "Up" | pi_N2hrsig == "Up" ~ "Up",  pi_N1hrsig == "NS" & pi_N2hrsig == "NS" ~ "NS")) %>%
  mutate(LY1 = case_when(pi_0.5hrsig == "Down" | pi_1hrsig == "Down" | pi_2hrsig == "Down" ~ "Down", pi_0.5hrsig == "Up" | pi_1hrsig == "Up" | pi_2hrsig == "Up" ~ "Up", pi_0.5hrsig == "NS"  & pi_1hrsig == "NS" & pi_2hrsig == "NS" ~ "NS")) %>%
  mutate(DHL4 = case_when(pi_D0.5hrsig == "Down" | pi_D1hrsig == "Down" | pi_D2hrsig == "Down" ~ "Down", pi_D0.5hrsig == "Up" | pi_D1hrsig == "Up" | pi_D2hrsig == "Up" ~ "Up", pi_D0.5hrsig == "NS"  & pi_D1hrsig == "NS" & pi_D2hrsig == "NS" ~ "NS"))  %>%
pivot_longer(NUD:DHL4, names_to = "cell_line", values_to = "group") %>%
  dplyr::group_by(cell_line, gbc_class, group) %>%
  dplyr::summarise(genes = n()) %>%
  dplyr::mutate(pos = paste(cell_line, gbc_class, sep = "_")) %>%
  dplyr::ungroup() %>%
  add_row(cell_line = "DHL4", gbc_class = "gbu", group = "Up", genes = 0, pos = "DHL4_gbu") %>%
  add_row(cell_line = "LY1", gbc_class = "gbd", group = "Down", genes = 0, pos = "LY1_gbd") %>%
  add_row(cell_line = "LY1", gbc_class = "gbu", group = "Up", genes = 0, pos = "LY1_gbu") %>%
  add_row(cell_line = "NUD", gbc_class = "gbu", group = "Up", genes = 0, pos = "NUD_gbu")

write_delim(siggenes_wpichange_cluster_count_0.5, "20230906 NUD_LY1_PS_pichanges_LY10.5hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_1, "20230906 NUD_LY1_PS_pichanges_LY11hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_2, "20230906 NUD_LY1_PS_pichanges_LY12hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_N1, "20230906 NUD_LY1_PS_pichanges_NUD1hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_N2, "20230906 NUD_LY1_PS_pichanges_NUD2hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_0.5, "20230906 NUD_LY1_PS_pichanges_DHL40.5hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_1, "20230906 NUD_LY1_PS_pichanges_DHL41hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_cluster_count_2, "20230906 NUD_LY1_PS_pichanges_DHL42hr_bycluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_nocluster, "20230906 NUD_LY1_PS_pichanges_nocluster_km3.txt", delim = "\t")
write_delim(siggenes_wpichange_gbc, "20230906 NUD_LY1_PS_pichanges_nocluster_gbc.txt", delim = "\t")
```

