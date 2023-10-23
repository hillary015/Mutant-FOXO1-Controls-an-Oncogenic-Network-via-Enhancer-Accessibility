# S6A. Annotation of regulated genes to proximal functional peak
The following code was used to annotate regulated genes to the closest functional peak in bash and R.

## Annotation
Regulated genes were annotated to the closest functional peak using HOMER in bash.
```
annotatePeaks.pl 20230622_ATACdown_eRNAdown_FKHbound_merge.txt hg19 -gtf 20230622_NUD_LY1_PSsigdown_refseq.gtf.filepart > 20230622_ATACdown_eRNAdown_FKHbound_merge_annotate_PSsiggenes_NUD_LY1_overlap.txt
```
## Bar graph
A bar graph show percentage of regulated genes with the closest functional peak within a given distance was generated in R.
### Load Libraries
```{r}
library(tidyverse)
library(hillaryscolors)
```

### Load files
```{r}
ATAC_down <- readr::read_delim("20230622_ATACdown_eRNAdown_FKHbound_merge_annotate_PSsiggenes_NUD_LY1_overlap.txt")
gene_IDs <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE//Comb_seqruns/20220809 LY1_FOXO1_PS_expressedgenes_duprm.txt")
```

### Process data for histogram
```{r}
options(scipen = 99)

funcfkh_bed <- ATAC_down %>%
  dplyr::filter(`Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt") %>%
select(Chr, Start, End) %>%
  write_delim("20230703_funcFKHmotifs.bed", delim = "\t")

ATAC_down_select <- ATAC_down %>%
  dplyr::filter(`Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt") %>%
  dplyr::mutate(Distance_to_nearest_TSS_kb = abs((as.numeric(`Distance to TSS`/1000)))) %>%
  dplyr::select(PeakID = `PeakID (cmd=annotatePeaks.pl 20230622_ATACdown_eRNAdown_FKHbound_merge.txt hg19 -gtf 20230622_NUD_LY1_PSsigdown_refseq.gtf.filepart)`, Distance_to_nearest_TSS_kb) %>%
  dplyr::mutate(bin = cut(Distance_to_nearest_TSS_kb, breaks = seq(0, 150000, 1), right = FALSE)) %>%
  tidyr::extract(bin, c("binstart", "binend"), "([[:alnum:]]+),([[:alnum:]]+)")


ATAC_down_select$binstart <- as.numeric(ATAC_down_select$binstart)
ATAC_down_select$binend <- as.numeric(ATAC_down_select$binend)
readr::write_delim(ATAC_down_select, "20230622_funcFKHmotifs_annotoPSgenes_histbins.txt", delim = "\t")

ATAC_down_bin_summ <- ATAC_down_select %>%
  dplyr::group_by(grp = cut(binstart, breaks = seq(0, 150000, 100))) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percent = (count/380)*100)
readr::write_delim(ATAC_down_bin_summ, "20230622_funcFKHmotifs_annotoPSgenes_histbins_summ.txt", delim = "\t")

ATAC_down_bin_summ_1MB <- ATAC_down_select %>%
  dplyr::group_by(grp = cut(binstart, c(0, 1000, -Inf))) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percent = (count/380)*100)
readr::write_delim(ATAC_down_bin_summ_1MB, "20230622_funcFKHmotifs_annotoPSgenes_histbins_summ_1MB.txt", delim = "\t")

ATAC_down_bin_summ_50 <- ATAC_down_select %>%
  dplyr::group_by(grp = cut(binstart, breaks = seq(0, 150000, 50))) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percent = (count/380)*100)
readr::write_delim(ATAC_down_bin_summ_50, "20230622_funcFKHmotifs_annotoPSgenes_histbins_50kb_summ.txt", delim = "\t")

ATAC_down_bin_summ_25 <- ATAC_down_select %>%
  dplyr::group_by(grp = cut(binstart, breaks = seq(0, 150000, 25))) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percent = (count/380)*100)
readr::write_delim(ATAC_down_bin_summ_25, "20230622_funcFKHmotifs_annotoPSgenes_histbins_25kb_summ.txt", delim = "\t")
```

### ATAC distance by gene
```{r}
Transcript <- c("NM_001289139", "NM_015184", "NM_001350898", "NM_033181", "NM_001308243", "NM_001364867", "NM_001242630", "NM_001142598", "NM_001385621", "NM_001328619", "NM_001394959", "NM_001322924", "NM_001377938", "NM_001243186", "NM_001354998", "NM_001371277", "NM_002399", "NM_001353625", "NM_033285", "NM_001329964", "NM_012418", "NM_001286761", "NM_003379", "NM_001079874", "NR_120485", "NM_001318327", "NM_001354817", "NM_001330491", "NM_001354642", "NM_145017", "NM_001385620", "NM_001350900", "NM_138810", "NM_001304481", "NM_001270765", "NR_033258", "NM_001278359", "NM_001291866", "NM_001370298", "NR_134961", "NM_001244638", "NM_001008540", "NM_001382691","NM_001350896", "NM_001257291", "NM_001354812", "NM_001286762", "NR_126060", "NR_146938", "NM_001348056", "NM_001376591", "NM_001370546", "NM_001348512", "NM_001308242", "NM_001321813", "NR_134960", "NM_001286688", "NM_001377936", "NM_001286691", "NM_001329145", "NM_001039570", "NM_001377345", "NM_007289", "NM_001377278", "NM_001258038", "NR_136519", "NM_001242629", "NM_001003688", "NM_005841", "NR_149144", "NM_001330374", "NM_004199", "NM_001377349", "NM_001193338", "NM_000626", "NM_001313915", "NM_021813", "NM_001149", "NM_001376256", "NM_001243999", "NM_001888", "NM_014795", "NM_001371275", "NR_026691", "NM_001369882", "NM_020987", "NM_001363551", "NM_001396408", "NM_001365677", "NM_001367282", "NM_001244", "NM_001278733", "NM_001382675", "NM_007333", "NM_005025", "NM_001297611", "NM_001330343", "NR_135810", "NM_001243786", "NM_001278368", "NM_001013255", "NM_001324222", "NM_001135196", "NM_001270764", "NR_104260", "NM_177532", "NM_001201404", "NM_173568", "NM_001319681", "NM_001270764", "NM_001135196", "NM_001198979", "NM_001278648", "NM_001162422", "NM_172315", "NM_001329804", "NM_001172477", "NM_001190259", "NM_001257407", "NM_001330663", "NM_014863", "NM_001329798", "NM_198196", "NM_016150", "NM_001174072", "NM_001145436", "NM_001396409", "NM_001328621", "NM_145341", "NM_001004416", "NM_001130084", "NM_198679", "NM_001278911", "NM_001323375", "NM_001331034", "NM_001204404", "NM_001195053", "NM_001198980", "NM_005080", "NM_001173982", "NM_001257997", "NM_000536"," NM_001017974", "NM_001135585", "NR_134959", "NM_001354997", "NM_001350901", "NR_134962", "NM_001270392", "NM_001282616", "NM_001017974", "NM_001286262", "NM_001377277", "NM_152785", "NM_001385619", "NM_001329150", "NR_138140", "NM_054114", "NR_073189")

Gene <- c("C14orf32", "PLCL2", "PION", "CNR1", "IL2RA", "IFI16", "GFOD1", "P4HA2", "ZNF608", "SLC2A5", "MARCH1", "RCSD1", "RAPGEF1", "PIM1", "SCL35E3", "IKZF2", "MEIS2", "C12orf49", "TP53INP1", "TP63", "FSCN2", "KIAA0226L", "EZR", "VAV3", "IL7R", "DNMBP", "SMAD1", "BNIP3L", "MME", "FLJ32771", "ZNF608", "PION", "TAGAP", "FGD4", "CHST15", "ZEB2", "RHOH", "GALNT2", "FGD4", "MYB", "ARID5B", "CXCR4", "MBNL2", "PION", "SCL9A7", "SMAD1", "KIAA0226L", "SERIN5C", "PION", "CXCR4", "IFI16", "CNR1", "CCDC85A", "IL2RA", "CBLB", "MYB", "ABLIM2", "RAPGEF1", "MLLT3", "TP63", "KREMEN1", "LRIG1", "MME", "RAG1", "SPRY1", "RCSD1", "GFOD1",  "SMAD1", "SPRY1", "SLC35E3", "FGD4", "P4HA2", "LRIG1", "FAIM3", "CD79B", "EIF2AK3", "BACH2", "ANK3", "CRYM", "SPIB", "CRYM", "ZEB2", "IKZF2", "IL10RA", "ADM2", "ANK3", "PLA2G15", "PLAC2", "P4HA2", "MEF2B", "TNFSF8", "TAGAP", "MBNL2", "KLRC3", "SERPINI1", "TMCC2","AICDA", "CBLB", "RAG2", "RHOH", "LSP1", "HMGCS1", "C10orf71", "CHST15", "PRDM15", "RASSF6", "WASF2", "UMODL1", "TCTN1", "CHST15", "C10orf71", "SMAP2", "ZBTB48", "ETS1", "MEIS2", "PLEKHG1", "RRM2B", "GCET2", "IL4R", "HMGCS1", "CHST15", "PLEKHG1", "CD96", "ASB2", "SERINC5", "LSS", "PLAC2", "SLC2A5", "PDCD4", "UMODL1", "ABLIM2", "RAPGEF2", "KCNMB2", "CAMK4", "TMCC2", "ANK3", "DDIT3", "SMAP2", "XBP1", "CHST11", "IL4R", "RAG2", "P4HA2", "SLC2A5", "MYB", "SLC35E3", "PION", "MYB", "RASSF6", "BEST3", "P4HA2", "TCP11L2", "RAG1", "GCET2", "ZNF608", "TP63", "PLEKHG1", "TAGAP", "HRK")

missannotated <- dplyr::data_frame(Transcript, Gene)
gene_IDs_adj <- gene_IDs %>%
  dplyr::add_row(missannotated)

ATAC_down_bygene <- ATAC_down %>%
    dplyr::filter(`Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt" | `Focus Ratio/Region Size` == "../FOXO1_eRNA_analysis/20230618_LY1_Nud_DHL4_eRNA_downpeaks_merged.txt|../Nud_DHL4_ATAC_forTobias/TF_footprinting_out/20230618_expressed_FKH_0hboundmtifs_merged_sorted_cat_bed4.txt") %>%
  dplyr::mutate(Distance_to_nearest_TSS_kb = abs((`Distance to TSS`/1000))) %>%
  dplyr::left_join(gene_IDs_adj, by = c("Entrez ID" = "Transcript")) %>%
  dplyr::select(PeakID = `PeakID (cmd=annotatePeaks.pl 20230622_ATACdown_eRNAdown_FKHbound_merge.txt hg19 -gtf 20230622_NUD_LY1_PSsigdown_refseq.gtf.filepart)`, Distance_to_nearest_TSS_kb, Gene) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(closest_ATAC_kb = min(Distance_to_nearest_TSS_kb))
readr::write_delim(ATAC_down_bygene, "20230623 func_FKHmotifs_losestATACdown_toPSgenedown.txt", delim = "\t")

ATAC_down_bygene_histogram <- ATAC_down_bygene %>%
  dplyr::mutate(bin = case_when(closest_ATAC_kb < 25 ~ "<25", closest_ATAC_kb < 50 ~ "<50", closest_ATAC_kb < 100 ~ "<100", closest_ATAC_kb < 250 ~ "<250", closest_ATAC_kb < 500 ~ "<500", closest_ATAC_kb < 1000 ~ "<1000", closest_ATAC_kb >= 1000 ~ ">=1000")) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percent = ((count/62) * 100)) %>%
  dplyr::ungroup() %>%
  dplyr::add_row(bin = "noATAC", count = 3, percent = ((3/62)*100))
readr::write_delim(ATAC_down_bygene_histogram, "20230623 func_FKHmotifs_closestATACdown_toPSgenedown_binsum.txt", delim = "\t")
```

### Histogram
```{r}
position <- c("noATAC", ">=1000", "<1000", "<500", "<250", "<100", "<50", "<25")
hist_bygene <- ggplot (ATAC_down_bygene_histogram, aes(x = bin, y = percent)) +
  geom_bar(stat ="identity", position = position_dodge(0.65), color = "black", fill = "#AFAFAF", width = 0.5) +
  ggtitle("Distance to nearest start site") +
  labs(x = "Distance from Peak Center", y = "Peaks per bin") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_discrete(limits = position) +
  theme_hillary() +
  coord_cartesian(clip = "off") +
  coord_flip() +
  theme(plot.title = element_text(vjust = 1, size = "20", face = "bold"))
hist_bygene
```

### Print to png
```{r}
Cairo::Cairo(file = "20230623 func_FKHmotifs_losestATACdowntoPSgene_4x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 4, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
hist_bygene
dev.off()
```
