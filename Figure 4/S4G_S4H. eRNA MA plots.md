# S4G_S4H eRNA MA plot
The following code is an example of the code used to make MA plots for intergenic RNA from PRO-seq data

## Load libraries
```{r}
library(tidyverse)
library(hillaryscolors)
library(viridis)
```

## Load data
```{r}
enhancer_change_1 <- readr::read_delim("20230307_Enhancer_change_T0vT1_deseq.txt")
enhancer_change_2 <- readr::read_delim("20230307_Enhancer_change_T0vT2_deseq.txt")

T0v1 <- readr::read_delim("count_enhancer_0v1.txt")
T0v2 <- readr::read_delim("count_enhancer_0v2.txt")
```

## Add threhold column
```{r}
thresholdfun <- function(data) {
   data %>%
    mutate(threshold = case_when(padj >= 0.05 ~ "NS", log2FoldChange > -0.585 & log2FoldChange < 0.585 ~ "NS", padj < 0.05 & log2FoldChange > 0.585 ~ "Up", padj < 0.05 & log2FoldChange < -0.585 ~ "Down"))
}

T1_thresh <- thresholdfun(enhancer_change_1)
T2_thresh <- thresholdfun(enhancer_change_2)

thresholdsumm <- function(data) {
  data %>%
    group_by(threshold) %>%
    summarise(count = n())
}

T1_thresh_summ <- thresholdsumm(T1_thresh)
print(T1_thresh_summ)
T2_thresh_summ <- thresholdsumm(T2_thresh)
print(T2_thresh_summ)
```

## Make plots
```{r}
MA_plot <- function(data, title, pal) {
  ggplot(data %>% arrange(desc(data)), aes(x = baseMean, y = log2FoldChange, colour = threshold)) +
  geom_point(size = 1) +
  ggtitle(title) +
  ylab("log2 fold change") +
  xlab("log10(baseMean)") +
  scale_y_continuous(breaks = seq(-4, 4, by = 2), expand = c(0, 0), limits = c(-4, 4)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 1, colour = "#434343") +
  scale_x_log10(limits = c(10, 100000)) +
  theme_hillary() +
  scale_color_manual(values = pal) +
  coord_cartesian(clip = "off") +
  annotation_logticks(sides = "b", outside = TRUE, size = 1) +
  theme(plot.title = element_text(vjust = 2)) 
}

updownpal <- c("#31607f", "#B3B6B7", "#70000C")

MA_1 <- MA_plot(T1_thresh, "FOXO1-FKBP 0v1", updownpal)
MA_1

MA_2 <- MA_plot(T2_thresh, "FOXO1-FKBP 0v2", updownpal)
MA_2
```

## Make peak files
```{r}
T1_wloc <- T0v1  %>%
    left_join(T1_thresh, by = c( "Enhancer_ID" = "row"))

T2_wloc <- T0v2  %>%
    left_join(T2_thresh, by = c( "Enhancer_ID"="row" ))

makepeak <- function(df, file, thr) {
df  %>%
  dplyr::filter(threshold == thr) %>%
  dplyr::select(Enhancer_ID, chr, start, end) %>%
  dplyr::mutate(strand = 0) %>%
  write.table(file = file, append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
}
makepeak(T1_wloc, "NUD_FOXO1_PS_eRNAdownpeak_1.txt", "Down")
makepeak(T2_wloc, "NUD_FOXO1_PS_eRNAdownpeak_2.txt", "Down")
makepeak(T1_wloc, "NUD_FOXO1_PS_eRNAuppeak_1.txt", "Up")
makepeak(T2_wloc, "NUD_FOXO1_PS_eRNAuppeak_2.txt", "Up")

makepeaknothr <- function(df, file) {
  df %>%
    dplyr::select(Enhancer_ID, chr, start, end) %>%
    dplyr::mutate(strand = 0) %>%
    write.table(file = file, append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
}

makepeaknothr(T1_wloc, "NUD_FOXO1_PS_eRNApeak_1.txt")
makepeaknothr(T2_wloc, "NUD_FOXO1_PS_eRNApeak_2.txt")
```

## Print to png
```{r}
cairopng_6x5 <- function(title, object) {
 Cairo::Cairo(file = title, 
             bg = "white",
             type = "png",
             units = "in", 
             width = 6, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
object
}

cairopng_6x5("20230314_NUDUL1_FOXO1_PS_eRNA_1_MAplots.png", MA_1)
dev.off()
cairopng_6x5("2020314_NUDUL1_FOXO1_PS_eRNA_2_MAplots.png", MA_2)
dev.off()
```


