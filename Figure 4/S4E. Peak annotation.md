# 4SE. Peak annotation
The following code is an example of the code used to annotate ATAC-seq peaks to the nearest genomic feature.

```
annotatePeaks.pl 20230822_LY1_NUD_DHL4_allpeaks_merge.txt hg19 > 20230822_LY1_NUD_DHL4_allpeaks_merge_annotate.txt

#        Peak file = 20230822_LY1_NUD_DHL4_allpeaks_merge.txt
#        Genome = hg19
#        Organism = human
#        Peak/BED file conversion summary:
#                BED/Header formatted lines: 0
#                peakfile formatted lines: 133925
#                Duplicated Peak IDs: 0

#        Peak File Statistics:
#                Total Peaks: 133925
#                Redundant Peak IDs: 0
#                Peaks lacking information: 0 (need at least 5 columns per peak)
#                Peaks with misformatted coordinates: 0 (should be integer)
#                Peaks with misformatted strand: 0 (should be either +/- or 0/1)

#        Peak file looks good!

#        Reading Positions...
#        -----------------------
#        Finding Closest TSS...
#        Annotating:........................
#                Annotation      Number of peaks Total size (bp) Log2 Ratio (obs/exp)    LogP enrichment (+values depleted)
#                3UTR    1913.0  26762632        0.724   -210.898
#                miRNA   16.0    96670   1.935   -11.677
#                TTS     2925.0  32286633        1.066   -645.937
#                pseudo  125.0   2089637 0.467   -8.004
#                Exon    3483.0  37006581        1.121   -841.978
#                Intron  66589.0 1258916630      0.289   -2228.981
#                Intergenic      47626.0 1691354378      -0.620  9876.128
#                Promoter        9795.0  35809872        2.660   -10080.314
#                5UTR    788.0   2586029 2.815   -867.554
#                snoRNA  0.0     357     -0.022  0.015
#                scRNA   0.0     97      -0.006  0.004
#        NOTE: If this part takes more than 2 minutes, there is a good chance
#                your machine ran out of memory: consider hitting ctrl+C and rerunning
#                the command with "-noann"
#        To capture annotation stats in a file, use "-annStats <filename>" next time
#        Annotating:........................
#                Annotation      Number of peaks Total size (bp) Log2 Ratio (obs/exp)    LogP enrichment (+values depleted)
#                3UTR    1913.0  26762632        0.724   -211.281
#                Other   21.0    3935577 -3.019  107.720
#                Unknown?        3.0     18134   1.935   -3.096
#                RNA     9.0     114374  0.863   -2.737
#                miRNA   16.0    96670   1.936   -11.684
#                ncRNA   665.0   6991458 1.136   -165.398
#                TTS     2925.0  32286633        1.066   -646.713
#                LINE    11450.0 622739191       -1.235  6771.107
#                LINE?   0.0     10448   -0.538  0.452
#                srpRNA  6.0     252957  -0.867  2.514
#                SINE    13710.0 378584817       -0.257  264.428
#                RC      15.0    442102  -0.351  1.577
#                tRNA    35.0    91526   3.144   -47.816
#                DNA?    7.0     264508  -0.709  2.148
#                pseudo  125.0   2089637 0.467   -8.023
#                DNA     3601.0  95640101        -0.200  40.546
#                Exon    3483.0  37006581        1.121   -842.935
#                Intron  42580.0 664015350       0.568   -3863.422
#                Intergenic      27670.0 863542451       -0.433  1857.951
#                Promoter        9793.0  35809872        2.660   -10080.687
#                5UTR    785.0   2586029 2.811   -862.026
#                snoRNA  0.0     357     -0.022  0.015
#                LTR?    0.0     21148   -0.937  0.915
#                scRNA   6.0     115467  0.264   -0.959
#                CpG-Island      2991.0  8588734 3.009   -3649.159
#                Low_complexity  245.0   15393271        -1.443  179.864
#                LTR     11170.0 259145993       -0.005  1.074
#                Simple_repeat   538.0   24806381        -0.996  168.161
#                snRNA   15.0    313162  0.147   -0.963
#                Unknown 70.0    1244230 0.379   -3.939
#                SINE?   3.0     42587   0.703   -1.270
#                Satellite       65.0    12339025        -3.038  335.677
#                rRNA    10.0    164138  0.494   -1.715
#        Counting Tags in Peaks from each directory...
#        Organism: human
#        Loading Gene Informaiton...
#        Outputing Annotation File...
#        Done annotating peaks file
```
## Make pie charts
### Load Libraries
```{r}
library(tidyverse)
library(viridis)
library(scales)
```

### Load annotated peaks file (from Homer) and extract simplified annotations to column anno
```{r}
ATAC_all_annotation <- read.delim("20230623_DHL4_LY1_Nud_allpeaks_merge_annotate.txt")

ATAC_all_annotation$anno = sapply(strsplit(x = as.character(ATAC_all_annotation$Annotation), split = ' (', fixed = T), '[[', 1)

ATAC_down_annotation <- read.delim("20230622_downpeaks_overlap_LY1_NUD_merge_bed_annotate.txt")

ATAC_down_annotation$anno = sapply(strsplit(x = as.character(ATAC_down_annotation$Annotation), split = ' (', fixed = T), '[[', 1)

ATAC_up_annotation <- read.delim("20230623_uppeaks_overlap_LY1_NUD_merge_bed_annotate.txt")

ATAC_up_annotation$anno = sapply(strsplit(x = as.character(ATAC_up_annotation$Annotation), split = ' (', fixed = T), '[[', 1)
```

### Summarize annotations
```{r}
anno_counts <- ATAC_all_annotation %>%
  group_by(anno) %>%
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))

anno_counts_down <- ATAC_down_annotation %>%
  group_by(anno) %>%
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))

anno_counts_up <- ATAC_up_annotation %>%
  group_by(anno) %>%
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))
```

### Write to txt file
```{r}
write.table(anno_counts, file = "20230622_annotatepeaks_ATACchange_all_annotated_summary.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)

write.table(anno_counts_down, file = "20230622_annotatepeaks_ATACchange_down_annotated_summary.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)

write.table(anno_counts, file = "20230622_annotatepeaks_ATACchange_all_annotated_summary.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)
```

### Make pie chart for anno
```{r}
pie <- ggplot(anno_counts, aes(x = "", y = n, fill = anno)) +
  ggtitle("ATAC all Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(8))
pie

piedown <- ggplot(anno_counts_down, aes(x = "", y = n, fill = anno)) +
  ggtitle("ATAC down Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(8))
piedown

pie_up <- ggplot(anno_counts_up, aes(x = "", y = n, fill = anno)) +
  ggtitle("ATAC UP Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(8))
pie_up
```

### Write to PNG
```{r}
Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_all_annotated_allanno.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
pie
dev.off()

Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_down_annotated_allanno.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
piedown
dev.off()

Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_up__annotated_allanno.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
pie_up
dev.off()
```

### Set all annotations besdies intergenic, intronic, and promoter to other 
```{r}
Other <- c("3' UTR", "5' UTR", "exon", "non-coding", "TTS")

anno_counts_other <- ATAC_all_annotation
anno_counts_other$anno2 <- sapply(anno_counts_other$anno, switch,
                  "3' UTR" = "Other", 
                  "5' UTR" = "Other", 
                  "exon" = "Other",
                  "Intergenic" = "Intergenic",
                  "intron" = "Intronic",
                  "non-coding" = "Other",
                  "promoter-TSS"	= "Promoter",
                  "TTS" = "Other")

anno_counts_other_sum <- anno_counts_other %>%
  group_by(anno2) %>%  
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))

anno_counts_other_down <- ATAC_down_annotation
anno_counts_other_down$anno2 <- sapply(anno_counts_other_down$anno, switch,
                  "3' UTR" = "Other", 
                  "5' UTR" = "Other", 
                  "exon" = "Other",
                  "Intergenic" = "Intergenic",
                  "intron" = "Intronic",
                  "non-coding" = "Other",
                  "promoter-TSS"	= "Promoter",
                  "TTS" = "Other")

anno_counts_other_sumdown <- anno_counts_other_down %>%
  group_by(anno2) %>%  
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))

anno_counts_other_up <- ATAC_up_annotation
anno_counts_other_up$anno2 <- sapply(anno_counts_other_up$anno, switch,
                  "3' UTR" = "Other", 
                  "5' UTR" = "Other", 
                  "exon" = "Other",
                  "Intergenic" = "Intergenic",
                  "intron" = "Intronic",
                  "non-coding" = "Other",
                  "promoter-TSS"	= "Promoter",
                  "TTS" = "Other")

anno_counts_other_sum_up <- anno_counts_other_up %>%
  group_by(anno2) %>%  
  summarise(n = n()) %>%
  mutate(percent = round(((n/sum(n))*100), digits = 2))
```

### Write to txt file
```{r}
write.table(anno_counts_other_sum, file = "20230622_annotatepeaks_ATACchange_all_annotated_summary_other.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)

write.table(anno_counts_other_sumdown, file = "20230622_annotatepeaks_ATACchange_down_annotated_summary_other.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)

write.table(anno_counts_other_sum_up, file = "20230622_annotatepeaks_ATACchange_up_annotated_summary_other.txt", sep = "\t",
            row.names = FALSE, quote =FALSE)
```

### Make pie chart for anno
```{r}
pie_other <- ggplot(anno_counts_other_sum, aes(x = "", y = n, fill = anno2)) +
  ggtitle("ATAC all Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(5))
pie_other

pie_otherdown <- ggplot(anno_counts_other_sumdown, aes(x = "", y = n, fill = anno2)) +
  ggtitle("ATAC down Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(5))
pie_otherdown

pie_other_up <- ggplot(anno_counts_other_sum_up, aes(x = "", y = n, fill = anno2)) +
  ggtitle("ATAC up Peak Annotations") + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = viridis(5))
pie_other_up
```

```{r}
Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_all_annotated_annoother.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
pie_other
dev.off()

Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_down_annotated_annoother.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
pie_otherdown
dev.off()

Cairo::Cairo(file="20230622_annotatepeaks_ATACchange_up_annotated_annoother.png", 
             bg="white",
             type="png",
             units="in", 
             width=5, 
             height=5, 
             pointsize=12, 
             dpi=300)
pie_other_up
dev.off()
```

