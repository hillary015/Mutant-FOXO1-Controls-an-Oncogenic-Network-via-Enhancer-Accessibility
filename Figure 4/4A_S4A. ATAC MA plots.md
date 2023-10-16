# 4A_S4A. ATAC MA Plots
The following code is an example of the code used to generate MA plots from ATAC-seq data.

##Load libraries
```{r}
library(tidyverse)
library(hillaryscolors)
```

##Load files
```{r}
ATACchange_0v15_deseqnorm <- read.delim("20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_15minv0hr.txt")
ATACchange_0v30_deseqnorm <- read.delim("20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_30v0hr.txt")
ATACchange_0v1_deseqnorm <- read.delim("20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_1hrv0hr.txt")
ATACchange_0v2_deseqnorm <- read.delim("20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_2hrv0hr.txt")
```

## Drop extra columns
```{r}
ATACchange_0v15_deseqnorm_dropna <- ATACchange_0v15_deseqnorm %>%
  select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
ATACchange_0v30_deseqnorm_dropna <- ATACchange_0v30_deseqnorm %>%
  select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
ATACchange_0v1_deseqnorm_dropna <- ATACchange_0v1_deseqnorm %>%
  select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
ATACchange_0v2_deseqnorm_dropna <- ATACchange_0v2_deseqnorm %>%
  select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
```

## Add threshhold column
```{r}
ATACchange_0v15_deseqnorm_dropna_thresh <- ATACchange_0v15_deseqnorm_dropna %>%
   mutate(thresh = case_when(padj < 0.05 & log2FoldChange <= -0.585 ~ "Down", pvalue >= 0.05 ~ "zNS", padj >= 0.05 ~ "zNS", log2FoldChange > -0.585 & log2FoldChange < 0.585 ~ "zNS", is.na(padj) ~ "zNS", padj < 0.05 & log2FoldChange >= -0.585 ~ "Up"))

ATACchange_0v30_deseqnorm_dropna_thresh <- ATACchange_0v30_deseqnorm_dropna %>%
  mutate(thresh = case_when(padj < 0.05 & log2FoldChange <= -0.585 ~ "Down", pvalue >= 0.05 ~ "zNS", padj >= 0.05 ~ "zNS", log2FoldChange > -0.585 & log2FoldChange < 0.585 ~ "zNS", is.na(padj) ~ "zNS", padj < 0.05 & log2FoldChange >= -0.585 ~ "Up"))

ATACchange_0v1_deseqnorm_dropna_thresh <- ATACchange_0v1_deseqnorm_dropna %>%
  mutate(thresh = case_when(padj < 0.05 & log2FoldChange <= -0.585 ~ "Down", pvalue >= 0.05 ~ "zNS", padj >= 0.05 ~ "zNS", log2FoldChange > -0.585 & log2FoldChange < 0.585 ~ "zNS", is.na(padj) ~ "zNS", padj < 0.05 & log2FoldChange >= -0.585 ~ "Up"))

ATACchange_0v2_deseqnorm_dropna_thresh <- ATACchange_0v2_deseqnorm_dropna %>%
  mutate(thresh = case_when(padj < 0.05 & log2FoldChange <= -0.585 ~ "Down", pvalue >= 0.05 ~ "zNS", padj >= 0.05 ~ "zNS", log2FoldChange > -0.585 & log2FoldChange < 0.585 ~ "zNS", is.na(padj) ~ "zNS", padj < 0.05 & log2FoldChange >= -0.585 ~ "Up"))
```

## Count each thresh
```{r}
thresh_counts_15_de <- ATACchange_0v15_deseqnorm_dropna_thresh %>%
  group_by(thresh) %>%
  summarise(n = n())
thresh_counts_15_de

thresh_counts_30_de <- ATACchange_0v30_deseqnorm_dropna_thresh %>%
  group_by(thresh) %>%
  summarise(n = n())
thresh_counts_30_de

thresh_counts_1_de <- ATACchange_0v1_deseqnorm_dropna_thresh %>%
  group_by(thresh) %>%
  summarise(n = n())
thresh_counts_1_de

thresh_counts_2_de <- ATACchange_0v2_deseqnorm_dropna_thresh %>%
  group_by(thresh) %>%
  summarise(n = n())
thresh_counts_2_de
```

## Make MA plots
```{r}
MA_0v15_de <- ggplot(ATACchange_0v15_deseqnorm_dropna_thresh %>% arrange(desc(thresh)), aes(x = baseMean, y = log2FoldChange, colour = thresh)) +
  geom_point(size = 0.5) +
  ggtitle("NUDUL1 FOXO1 ATAC 15 min") +
  ylab("log2 fold change") +
  xlab("log10(BaseMean)") +
  scale_y_continuous(breaks = seq(-6, 6, by = 3), expand = c(0, 0), limits = c(-6, 6)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1, colour = "#636363") +
  scale_x_log10(limits = c(1, 100000)) +
  theme_hillary() +
  scale_color_manual(values = c("#31607f", "#70000C", "#B3B6B7")) +
  annotation_logticks(sides = "b", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2)) 
MA_0v15_de

MA_0v30_de <- ggplot(ATACchange_0v30_deseqnorm_dropna_thresh %>% arrange(desc(thresh)), aes(x = baseMean, y = log2FoldChange, colour = thresh)) +
  geom_point(size = 0.5) +
  ggtitle("NUDUL1 FOXO1 ATAC 30 min") +
  ylab("log2 fold change") +
  xlab("log10(BaseMean)") +
  scale_y_continuous(breaks = seq(-6, 6, by = 3), expand = c(0, 0), limits = c(-6, 6)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1, colour = "#636363") +
  scale_x_log10(limits = c(1, 100000)) +
  theme_hillary() +
  scale_color_manual(values = c("#31607f", "#70000C", "#B3B6B7")) +
  annotation_logticks(sides = "b", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2)) 
MA_0v30_de

MA_0v1_de <- ggplot(ATACchange_0v1_deseqnorm_dropna_thresh %>% arrange(desc(thresh)), aes(x = baseMean, y = log2FoldChange, colour = thresh)) +
  geom_point(size = 0.5) +
  ggtitle("NUDUL1 FOXO1 ATAC 1 hr") +
  ylab("log2 fold change") +
  xlab("log10(BaseMean)") +
  scale_y_continuous(breaks = seq(-6, 6, by = 3), expand = c(0, 0), limits = c(-6, 6)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1, colour = "#636363") +
  scale_x_log10(limits = c(1, 100000)) +
  theme_hillary() +
  scale_color_manual(values = c("#31607f", "#70000C", "#B3B6B7")) +
  annotation_logticks(sides = "b", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2))  
MA_0v1_de

MA_0v2_de <- ggplot(ATACchange_0v2_deseqnorm_dropna_thresh %>% arrange(desc(thresh)), aes(x = baseMean, y = log2FoldChange, colour = thresh)) +
  geom_point(size = 0.5) +
  ggtitle("NUDUL1 FOXO1 ATAC 2 hr") +
  ylab("log2 fold change") +
  xlab("log10(BaseMean)") +
  scale_y_continuous(breaks = seq(-6, 6, by = 3), expand = c(0, 0), limits = c(-6, 6)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1, colour = "#636363") +
  scale_x_log10(limits = c(1, 100000)) +
  theme_hillary() +
  scale_color_manual(values = c("#31607f", "#70000C", "#B3B6B7")) +
  annotation_logticks(sides = "b", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2))  
MA_0v2_de
```

## Saves as png
```{r}
Cairo::Cairo(file ="20230509_NUDUL1_FOXO1_ATAC_deseqnorm_hg19chrMrm_0v15min_genrich_MA_L2FC_padjco_6x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 6, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
MA_0v15_de
dev.off()

Cairo::Cairo(file ="20230509_NUDUL1_FOXO1_ATAC_deseqnorm_hg19chrMrm_0v30min_genrich_MA_L2FC_padjco_6x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 6, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
MA_0v30_de
dev.off()

Cairo::Cairo(file ="20230509_NUDUL1_FOXO1_ATAC_deseqnorm_hg19chrMrm_0v1hr_genrich_MA_L2FC_padjco_6x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 6, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
MA_0v1_de
dev.off()


Cairo::Cairo(file ="20230509_NUDUL1_FOXO1_ATAC_deseqnorm_hg19chrMrm_0v2hr_genrich_MA_L2FC_padjco_6x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 6, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
MA_0v2_de
dev.off()
```



