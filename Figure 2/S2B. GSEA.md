# S2B. GSEA
This code was used to perform GSEA and generate a heatmap for PRO-seq data in R.

## Project description:
#### Use clusterprofiler to ifentify enriched gene signatures 

## Load Library
```{r}
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(hillaryscolors)
library(xlsx)
library(DOSE)
library(viridis)
```

## Load Data
```{r}
allexpressed_NUD<- readr::read_delim("20230205 NUDUL1_FOXO1_PS_expressedgenes_duprm_manfilter.txt")

NUD_0v1_rnk <- readr::read_delim("gb_change_0v1h.txt") %>%
  dplyr::filter(Transcript %in% allexpressed_NUD$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

NUD_0v2_rnk <- readr::read_delim("gb_change_0v2.txt") %>%
  dplyr::filter(Transcript %in% allexpressed_NUD$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

allexpressed_LY1 <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/20220912 LY1_FOXO1_PS_expressedgenes_duprm_manfilter.txt")

LY1_0v0.5_rnk <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/gb_change_0v0.5.txt") %>%
  dplyr::filter(Transcript %in% allexpressed_LY1$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

LY1_0v1_rnk <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/gb_change_0v1.txt") %>%
  dplyr::filter(Transcript %in% allexpressed_LY1$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

LY1_0v2_rnk <- readr::read_delim("../Exp 0358  PS_LY1_FOXO1_JE/Comb_seqruns/gb_change0v2.txt") %>%
  dplyr::filter(Transcript %in% allexpressed_LY1$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

allexpressed_DHL4 <- readr::read_delim("../Exp 0071/20230530 DHL4_FOXO1_PS_expressedgenes_duprm_manfilter.txt")

DHL4_0v0.5_rnk <- readr::read_csv("../Exp 0071/DHL4_T0vT0.5_gb_change_2.csv") %>%
  dplyr::filter(Transcript %in% allexpressed_DHL4$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

DHL4_0v1_rnk <- readr::read_csv("../Exp 0071/DHL4_T0vT1_gb_change_2.csv") %>%
  dplyr::filter(Transcript %in% allexpressed_DHL4$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

DHL4_0v2_rnk <- readr::read_csv("../Exp 0071/DHL4_T0vT2_gb_change_2.csv") %>%
  dplyr::filter(Transcript %in% allexpressed_DHL4$Transcript) %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(Gene, log2FoldChange)

victora <- readr::read_delim("../Exp 0131_FOXO1_ETO_RS/Victora_2012_LZ_DZ_gs.txt")
roberto <- readr::read_delim("Roberto_genesets.txt")
```

## Make gene sets
```{r}
#Use mSigDB hallmarks
hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(hallmarks)
#Custom Bcell sets
victora$gczone <- paste("Victora", victora$`Enriched in`, sep = "_")
victora_adj <- victora %>%
  dplyr::select(gczone, `Gene Symbol`) %>%
  dplyr::filter(`Gene Symbol` != "---") %>%
  dplyr::distinct()

roberto_adj <- roberto %>%
  dplyr::mutate(gene_set = case_when(...3 < 0 ~ "Down_in_mut", ...3 > 0 ~ "Up_in_mut")) %>%
  dplyr::select(gene_set, Gene = 1)
```

## Make ranked lists
```{r}
NUD_0v1_rnk_adj <- NUD_0v1_rnk$log2FoldChange
names(NUD_0v1_rnk_adj) <- NUD_0v1_rnk$Gene
genelist_NUD_1 = sort(NUD_0v1_rnk_adj, decreasing = TRUE)

NUD_0v2_rnk_adj <- NUD_0v2_rnk$log2FoldChange
names(NUD_0v2_rnk_adj) <- NUD_0v2_rnk$Gene
genelist_NUD_2 = sort(NUD_0v2_rnk_adj, decreasing = TRUE)

LY1_0v0.5_rnk_adj <- LY1_0v0.5_rnk$log2FoldChange
names(LY1_0v0.5_rnk_adj) <- LY1_0v0.5_rnk$Gene
genelist_LY1_0.5 = sort(LY1_0v0.5_rnk_adj, decreasing = TRUE)

LY1_0v1_rnk_adj <- LY1_0v1_rnk$log2FoldChange
names(LY1_0v1_rnk_adj) <- LY1_0v1_rnk$Gene
genelist_LY1_1 = sort(LY1_0v1_rnk_adj, decreasing = TRUE)

LY1_0v2_rnk_adj <- LY1_0v2_rnk$log2FoldChange
names(LY1_0v2_rnk_adj) <- LY1_0v2_rnk$Gene
genelist_LY1_2 = sort(LY1_0v2_rnk_adj, decreasing = TRUE)

DHL4_0v0.5_rnk_adj <- DHL4_0v0.5_rnk$log2FoldChange
names(DHL4_0v0.5_rnk_adj) <- DHL4_0v0.5_rnk$Gene
genelist_DHL4_0.5 = sort(DHL4_0v0.5_rnk_adj, decreasing = TRUE)

DHL4_0v1_rnk_adj <- DHL4_0v1_rnk$log2FoldChange
names(DHL4_0v1_rnk_adj) <- DHL4_0v1_rnk$Gene
genelist_DHL4_1 = sort(DHL4_0v1_rnk_adj, decreasing = TRUE)

DHL4_0v2_rnk_adj <- DHL4_0v2_rnk$log2FoldChange
names(DHL4_0v2_rnk_adj) <- DHL4_0v2_rnk$Gene
genelist_DHL4_2 = sort(DHL4_0v2_rnk_adj, decreasing = TRUE)
```

## Gene set enrichmenet analysis
```{r}
gsea_hallmarks_NUD_0v1 <- GSEA(genelist_NUD_1, TERM2GENE = hallmarks, pvalueCutoff = 1)
gsea_hallmarks_NUD_0v2 <- GSEA(genelist_NUD_2, TERM2GENE = hallmarks, pvalueCutoff = 1)

gsea_victora_NUD_0v1 <- GSEA(genelist_NUD_1, TERM2GENE = victora_adj, pvalueCutoff = 1)
gsea_victora_NUD_0v2 <- GSEA(genelist_NUD_2, TERM2GENE = victora_adj, pvalueCutoff = 1)

gsea_roberto_NUD_0v1 <- GSEA(genelist_NUD_1, TERM2GENE = roberto_adj, pvalueCutoff = 1)
gsea_roberto_NUD_0v2 <- GSEA(genelist_NUD_2, TERM2GENE = roberto_adj, pvalueCutoff = 1)

readr::write_delim(as.data.frame(gsea_hallmarks_NUD_0v1), "20230907 NUD_FOXO1_PS_combseqruns_hallmarksgsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_hallmarks_NUD_0v2), "20230907 NUD_FOXO1_PS_combseqruns_hallmarksgsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_NUD_0v1), "20230907 NUD_FOXO1_PS_combseqruns_victoragsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_NUD_0v2), "20230907 NUD_FOXO1_PS_combseqruns_victoragsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_NUD_0v1), "20230907 NUD_FOXO1_PS_combseqruns_robertogsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_NUD_0v2), "20230907 NUD_FOXO1_PS_combseqruns_robertogsea_0v2.txt", delim = "\t")

xlsx::write.xlsx(gsea_hallmarks_NUD_0v1, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = FALSE, sheet = "0v1 hm")
xlsx::write.xlsx(gsea_hallmarks_NUD_0v2, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 hm")
xlsx::write.xlsx(gsea_victora_NUD_0v1, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 v")
xlsx::write.xlsx(gsea_victora_NUD_0v2, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 v")
xlsx::write.xlsx(gsea_roberto_NUD_0v1, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 r")
xlsx::write.xlsx(gsea_roberto_NUD_0v2, "20230907 NUD_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 r")

gsea_hallmarks_LY1_0v0.5 <- GSEA(genelist_LY1_0.5, TERM2GENE = hallmarks, pvalueCutoff = 1)
gsea_hallmarks_LY1_0v1 <- GSEA(genelist_LY1_1, TERM2GENE = hallmarks, pvalueCutoff = 1)
gsea_hallmarks_LY1_0v2 <- GSEA(genelist_LY1_2, TERM2GENE = hallmarks, pvalueCutoff = 1)

gsea_victora_LY1_0v0.5 <- GSEA(genelist_LY1_0.5, TERM2GENE = victora_adj, pvalueCutoff = 1)
gsea_victora_LY1_0v1 <- GSEA(genelist_LY1_1, TERM2GENE = victora_adj, pvalueCutoff = 1)
gsea_victora_LY1_0v2 <- GSEA(genelist_LY1_2, TERM2GENE = victora_adj, pvalueCutoff = 1)

gsea_roberto_LY1_0v0.5 <- GSEA(genelist_LY1_0.5, TERM2GENE = roberto_adj, pvalueCutoff = 1)
gsea_roberto_LY1_0v1 <- GSEA(genelist_LY1_1, TERM2GENE = roberto_adj, pvalueCutoff = 1)
gsea_roberto_LY1_0v2 <- GSEA(genelist_LY1_2, TERM2GENE = roberto_adj, pvalueCutoff = 1)

readr::write_delim(as.data.frame(gsea_hallmarks_LY1_0v0.5), "20230907 LY1_FOXO1_PS_combseqruns_hallmarksgsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_hallmarks_LY1_0v1), "20230907 LY1_FOXO1_PS_combseqruns_hallmarksgsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_hallmarks_LY1_0v2), "20230907 LY1_FOXO1_PS_combseqruns_hallmarksgsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_LY1_0v0.5), "20230907 LY1_FOXO1_PS_combseqruns_victoragsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_LY1_0v1), "20230907 LY1_FOXO1_PS_combseqruns_victoragsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_LY1_0v2), "20230907 LY1_FOXO1_PS_combseqruns_victoragsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_LY1_0v0.5), "20230907 LY1_FOXO1_PS_combseqruns_robertogsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_LY1_0v1), "20230907 LY1_FOXO1_PS_combseqruns_robertogsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_LY1_0v2), "20230907 LY1_FOXO1_PS_combseqruns_robertogsea_0v2.txt", delim = "\t")

xlsx::write.xlsx(gsea_hallmarks_LY1_0v0.5, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = FALSE, sheet = "0v0.5 hm")
xlsx::write.xlsx(gsea_hallmarks_LY1_0v1, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 hm")
xlsx::write.xlsx(gsea_hallmarks_LY1_0v2, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 hm")
xlsx::write.xlsx(gsea_victora_LY1_0v0.5, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v0.5 v")
xlsx::write.xlsx(gsea_victora_LY1_0v1, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 v")
xlsx::write.xlsx(gsea_victora_LY1_0v2, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 v")
xlsx::write.xlsx(gsea_roberto_LY1_0v0.5, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v0.5 r")
xlsx::write.xlsx(gsea_roberto_LY1_0v1, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 r")
xlsx::write.xlsx(gsea_roberto_LY1_0v2, "20230907 LY1_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 r")

gsea_hallmarks_DHL4_0v0.5 <- GSEA(genelist_DHL4_0.5, TERM2GENE = hallmarks, pvalueCutoff = 1)
gsea_hallmarks_DHL4_0v1 <- GSEA(genelist_DHL4_1, TERM2GENE = hallmarks, pvalueCutoff = 1)
gsea_hallmarks_DHL4_0v2 <- GSEA(genelist_DHL4_2, TERM2GENE = hallmarks, pvalueCutoff = 1)

gsea_victora_DHL4_0v0.5 <- GSEA(genelist_DHL4_0.5, TERM2GENE = victora_adj, pvalueCutoff = 1)
gsea_victora_DHL4_0v1 <- GSEA(genelist_DHL4_1, TERM2GENE = victora_adj, pvalueCutoff = 1)
gsea_victora_DHL4_0v2 <- GSEA(genelist_DHL4_2, TERM2GENE = victora_adj, pvalueCutoff = 1)

gsea_roberto_DHL4_0v0.5 <- GSEA(genelist_DHL4_0.5, TERM2GENE = roberto_adj, pvalueCutoff = 1)
gsea_roberto_DHL4_0v1 <- GSEA(genelist_DHL4_1, TERM2GENE = roberto_adj, pvalueCutoff = 1)
gsea_roberto_DHL4_0v2 <- GSEA(genelist_DHL4_2, TERM2GENE = roberto_adj, pvalueCutoff = 1)

readr::write_delim(as.data.frame(gsea_hallmarks_DHL4_0v0.5), "20230907 DHL4_FOXO1_PS_combseqruns_hallmarksgsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_hallmarks_DHL4_0v1), "20230907 DHL4_FOXO1_PS_combseqruns_hallmarksgsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_hallmarks_DHL4_0v2), "20230907 DHL4_FOXO1_PS_combseqruns_hallmarksgsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_DHL4_0v0.5), "20230907 DHL4_FOXO1_PS_combseqruns_victoragsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_DHL4_0v1), "20230907 DHL4_FOXO1_PS_combseqruns_victoragsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_victora_DHL4_0v2), "20230907 DHL4_FOXO1_PS_combseqruns_victoragsea_0v2.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_DHL4_0v0.5), "20230907 DHL4_FOXO1_PS_combseqruns_robertogsea_0v0.5.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_DHL4_0v1), "20230907 DHL4_FOXO1_PS_combseqruns_robertogsea_0v1.txt", delim = "\t")
readr::write_delim(as.data.frame(gsea_roberto_DHL4_0v2), "20230907 DHL4_FOXO1_PS_combseqruns_robertogsea_0v2.txt", delim = "\t")

xlsx::write.xlsx(gsea_hallmarks_DHL4_0v0.5, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = FALSE, sheet = "0v0.5 hm")
xlsx::write.xlsx(gsea_hallmarks_DHL4_0v1, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 hm")
xlsx::write.xlsx(gsea_hallmarks_DHL4_0v2, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 hm")
xlsx::write.xlsx(gsea_victora_DHL4_0v0.5, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v0.5 v")
xlsx::write.xlsx(gsea_victora_DHL4_0v1, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 v")
xlsx::write.xlsx(gsea_victora_DHL4_0v2, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 v")
xlsx::write.xlsx(gsea_roberto_DHL4_0v0.5, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v0.5 r")
xlsx::write.xlsx(gsea_roberto_DHL4_0v1, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v1 r")
xlsx::write.xlsx(gsea_roberto_DHL4_0v2, "20230907 DHL4_FOXO1_PS_combseqruns_Hallmarks_Victora_Roberto_clusterprofiler.xlsx", append = TRUE, sheet = "0v2 r")
```

