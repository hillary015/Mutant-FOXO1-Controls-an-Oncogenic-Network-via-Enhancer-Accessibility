# S1B. Pearson Correlation
The following code was used to calculate RNA-seq Pearson correlation values and graphs in R. 
## Load libraries
```{r}
library(tidyverse)
library(hillaryscolors)
```

## Load files
```{r}
NUD <- read_delim("NUD_parvfkbp_gene_exp_diff.txt")
LY7 <- read_delim("LY7_parvfkbp_gene_exp_diff.txt")
```


## Process file
```{r}
NUD_adj <- NUD %>%
  dplyr::filter(status == "OK") %>%
  select(gene_id, par = value_1, T0 = value_2)

LY7_adj <- LY7 %>%
  dplyr::filter(status == "OK") %>%
  select(gene_id, par = value_1, T0 = value_2)
```

## Correlation
```{r}
NUD_cor <- cor.test(NUD_adj$par, NUD_adj$T0, 
                    method = "pearson")
NUD_cor

LY7_cor <- cor.test(LY7_adj$par, LY7_adj$T0, 
                    method = "pearson")
LY7_cor
```

## Plot Graph
```{r}
scatter_NUD <- ggplot(NUD_adj, aes(x = log10(par), y = log10(T0))) +
  geom_point(size = 0.25, color = "#979797") +
  ggtitle("NUD All Expressed Genes (FKBP vs. WT)") +
  ylab("log10(taggedFPKM)") +
  xlab("log10(parentalFPKM)") +
  theme_hillary() +
  scale_x_continuous(expand = c(0, 0), limits = c(-2, 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 5)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2)) 
scatter_NUD

scatter_LY7 <- ggplot(LY7_adj, aes(x = log10(par), y = log10(T0))) +
  geom_point(size = 0.25, color = "#979797") +
  ggtitle("LY7 All Expressed Genes (FKBP vs. WT)") +
  ylab("log10(taggedFPKM)") +
  xlab("log10(parentalFPKM)") +
  theme_hillary() +
  scale_x_continuous(expand = c(0, 0), limits = c(-2, 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 5)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(vjust = 2)) 
scatter_LY7
```

## Print to png
```{r}
Cairo::Cairo(file = "20230509_NUD_FOXO1_RS_parvFKBP_FPKMcor.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
scatter_NUD
dev.off()

Cairo::Cairo(file = "20230509_LY7_FOXO1_RS_parvFKBP_FPKMcor.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
scatter_LY7
dev.off()
```
