# ATAC-seq Data Processing
The code used to process ATAC-seq data is presented below.
## Adaptor Trimming
Adaptor trimming was performed using Trimmomatic in bash. Please note that file names may need to be changed based on naming conventions at the time of sequencing.
```
for file in $(ls *_R1_001.fastq.gz); #Lists all R1.fastq files in directory
do
base=${file//_R1_001.fastq.gz} #Extracts base filename, ie XXXX_HL-X_S1_L005 
if [ ! -f "${base}"-forward_paired.fastq.gz ]; #Executes the following code if there is no corresponding index file
then
echo "${base}" #Prints base
java -classpath /Users/hiebertlab/software/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 "${base}"_R1_001.fastq.gz "${base}"_R2_001.fastq.gz "${base}"-forward_paired.fastq.gz "${base}"-forward_unpaired.fastq.gz "${base}"-reverse_paired.fastq.gz "${base}"-reverse_unpaired.fastq.gz ILLUMINACLIP:../activemotif_ATAC_adaptors.txt:2:30:7 LEADING:15 TRAILING:15 MINLEN:15 #Trims files
echo "Trimming of "${base}"_R1_001.fastq.gz "${base}"_R2_001.fastq.gz complete" #Prints string to console when trimming is complete
else
echo "Skipping "${base}" as "${base}"-forward_paired.fastq.gz already exists" #Prints string to console if corresponding BAM file already exists
fi
done
```
## Alignment
Trimmed data was aligned to a concatenated human (hg19)/ Drosophila (dm3) genome using Bowtie2 in bash.
```
for file in $(ls *-forward_paired.fastq.gz); #Lists all R1.fastq files in directory
do
if [ ! -f ${file//-forward_paired.fastq.gz}.hg19dm3.sam ]; #Executes the following code if SAM file doesn't exist
then
base=${file//-forward_paired.fastq.gz} #Extracts base filename, ie XXXX_HL-X_S1_L005
echo "${base}" #Prints base
bowtie2 -p 16 -x /Volumes/Hillary_3/genomes/hg19_dm3/hg19+dm3 -1 "${base}"-forward_paired.fastq.gz -2 "${base}"-reverse_paired.fastq.gz --no-mixed --no-discordant --phred33 -X 2000 -S "${base}".hg19dm3.sam 
echo "Alignment of -1 "${base}"-forward_paired.fastq.gz -2 "${base}"-reverse_paired.fastq.gz to hg19dm3 complete" #Prints string to console when trimming is complete
else
echo "Skipping alignment of "${file//-forward_paired.fastq.gz}" as "${file//-forward_paired.fastq.gz}".hg19dm3.sam already exists" #Prints string to console if SAM file already exists
fi
done

```
## SAM/BAM conversion and read filtering
Aligned SAM files were converted to BAM files and low quality reads were removed using SAMtools in bash.
```
for file in $(ls *.sam); #Creates a list of SAM files in working directory
do
if [ ! -f ${file%.*}.bam ]; #Executes following code if SAM does not have corresponding BAM
then
echo "Converting $file to ${file%.*}.bam" #Prints string to console
samtools view -S -b -@8 $file > ${file%.*}.bam # Converts SAM file to BAM file
else
echo "Skipping conversion of $file as ${file%.*}.bam already exists" #Prints string to console if corresponding BAM file already exists
fi

if [ ! -f ${file%.*}.sorted.bam ]; #Executes the following code if BAM file hasn't been sorted
then
echo "Sorting ${file%.*}.bam and creating new file ${file%.*}.sorted.bam" #Prints string to console
samtools sort -@8 -T temp ${file%.*}.bam -o ${file%.*}.sorted.bam  #Resorts BAM file
else
echo "Skipping sorting of ${file%.*}.bam as ${file%.*}.sorted.bam already exists" #Prints string to console if sorted BAM file already exists
fi
if [ ! -f ${file%.*}.F4q10.sorted.bam ]; #Executes the following code if BAM file hasn't been filtered
then
echo "Filtering low quality reads (F4, q10) from ${file%.*}.sorted.bam and creating new file ${file%.*}.F4q10.sorted.bam" #Prints to console
samtools view -b -F 4 -q 10 ${file%.*}.sorted.bam > ${file%.*}.F4q10.sorted.bam #Filters out low quality reads
else
echo "Skipping filtering of ${file%.*}.sorted.bam as ${file%.*}.F4q10.sorted.bam already exists" #Prints string to console if filtered BAM file already exists
fi

if [ ! -f ${file%.*}.F4q10.sorted.bam.bai ]; #Executes the following code if BAM file hasn't been indexed
then
echo "Indexing ${file%.*}.F4q10.sorted.bam and creating new file ${file%.*}.F4q10.sorted.bam.bai" #Prints string to console
samtools index ${file%.*}.F4q10.sorted.bam #Indexes BAM file
else
echo "Skipping indexing of ${file%.*}.F4q10.sorted.bam as ${file%.*}.F4q10.sorted.bam.bai already exists" #Prints string to console if indexed BAM file already exists
fi
done
```
## Count and remove spike in reads
Spike-in reads were counted and removed using SAMtools in bash. Mitochondrial reads were also removed.
```
for file in $(ls *.F4q10.sorted.bam); #Creates a list of sorted BAM files in working directory
do
if [ ! -f idxstats_${file%.*}.txt ]; #Executes following code if sorted.BAM does not have corresponding txt file
then
echo "Reporting alignment statistics for $file and printing to idxstats_${file%.*}.txt" #Prints string to console
samtools idxstats $file > idxstats_${file%.*}.txt #Reports alignment statistics for sorted.bam file and prints to .txt file
else
echo "Skipping $file as ${file%.*}.sorted.txt already exists" #Prints string to console if corresponding BAM file already exists
fi
done

## Make new file with only hg19 reads (-chrM)
for file in $(ls *.F4q10.sorted.bam); #Creates a list of sorted BAM files in working directory 
do
if [ ! -f ${file%.*}.chrMrm.bam ]; #Executes following code if mito reads have not yet been removed from sorted.bam file
then
echo "Printing hg19 reads from $file to ${file%.*}.chrMrm.bam" #Prints string to console
samtools view -bh $file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file%.*}.chrMrm.bam
else
echo "Skipping $file as ${file%.*}.chrMrm.bam already exists" #Prints string to console if corresponding BAM file already exists
fi
done

## Sorts new file
for file in $(ls *.chrMrm.bam); #Creates a list of BAM files with spikein reads removed (dm3rm) in working directory
do
if [ ! -f ${file%.*}.sorted.bam ]; #Executes following code if BAM file hasn't been resorted
then
echo "Sorting ${file%.*}.bam and creating new file ${file%.*}.sorted.bam" #Prints string to console
samtools sort -@ 10 -T ${file%.*}.temp -o ${file%.*}.sorted.bam ${file%.*}.bam #Reorts BAM file
else
echo "Skipping sorting of ${file%.*}.bam as ${file%.*}.sorted.bam already exists" #Prints string to console if resorted BAM file already exists
fi

## Indexes sorted file
if [ ! -f ${file%.*}.sorted.bam.bai ]; #Executes following code if BAM file hasn't been indexed
then
echo "Indexing ${file%.*}.sorted.bam and creating new file ${file%.*}.sorted.bam.bai" #Prints string to console
samtools index ${file%.*}.sorted.bam #Indexes BAM file
else
echo "Skipping indexing of ${file%.*}.sorted.bam as ${file%.*}.sorted.bam.bai already exists" #Prints string to console if indexed BAM file already exists
fi
done
```
## Name sort files for Genrich
Files were name-sorted for use in Genrich using SAMtools in bash.
```
for file in $(ls *.chrMrm.bam); #Creates a list of BAM files with spikein reads removed (dm3rm) in working directory
do
if [ ! -f ${file%.*}.namesorted.bam ]; #Executes following code if BAM file hasn't been resorted
then
echo "Sorting ${file%.*}.bam and creating new file ${file%.*}.namesorted.bam" #Prints string to console
samtools sort -@ 10 -n -T ${file%.*}.temp -o ${file%.*}.namesorted.bam ${file%.*}.bam #Resorts BAM file
else
echo "Skipping sorting of ${file%.*}.bam as ${file%.*}.namesorted.bam already exists" #Prints string to console if resorted BAM file already exists
fi
done
```

## Call peaks with Genrich
Peaks were called in using Genrich in conda.
```
Genrich -t /Volumes/Hillary_3/9330-subset/9330-HL-0001_S1_L005.hg19dm3.F4q10.sorted.chrMrm.namesorted.bam,/Volumes/Hillary_3/9330-subset/9330-HL-0002_S1_L005.hg19dm3.F4q10.sorted.chrMrm.namesorted.bam -j -r -E /Volumes/Hillary_3/genomes/hg19_blacklist/hg19.blacklistpeaks.bed -q 0.05 -o 9330-0hr_genrich.j.r.Eblacklist.q0.05.narrowPeak
```
## Quantification of reads and differential expression analysis
A consensus peakset was defined using DiffBind in R. Reads within the peaks defined in the consensus peakset were quantified with DiffBind in R. Counts were normalized and differential expression analysis was performed with DESeq2 in R.
```
## Project Description:
Use diffbind and DEseq2 to identify changes in chromatin accessibility loci after dTAG47 treatment

## Load libraries 
library("DiffBind")
packageVersion("DiffBind")
library("DESeq2")
packageVersion("DESeq2")
library(tidyverse)
library(gridExtra)
library(viridis)

## Load Data
##Data loaded from external hard drive (.bai files corresponding to .bam files loaded need to be located in the same directory), loaded into samplesheet format
ATAC_allnarrowpeaks <- data.frame(SampleID = c("0hrA", "0hrB", "15minA", "15minB", "30minA", "30minB", "1hrA", "1hrB", "2hrA", "2hrB"),
                                  Condition = c("0hr", "0hr", "15min", "15min", "30min", "30min", "1hr", "1hr", "2hr", "2hr"),
                                  Replicate = c("A", "B", "A", "B", "A", "B", "A", "B", "A", "B"),
                                  bamReads = c("E:/9330-NUDUL1_ATAC/9330-HL-0001_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0002_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0003_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0004_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0005_S1_L005_samp81mill.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0006_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0007_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0008_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0009_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam", "E:/9330-NUDUL1_ATAC/9330-HL-0010_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam"),
                                  Peaks = c("E:/9330-NUDUL1_ATAC/Genrich/9330-0hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-0hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-0.25hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-0.25hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-0.5hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-0.5hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-1hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-1hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-2hr_genrich.j.r.Eblacklist.q0.05.narrowPeak", "E:/9330-NUDUL1_ATAC/Genrich/9330-2hr_genrich.j.r.Eblacklist.q0.05.narrowPeak"), PeakCaller = "narrow")

## Diffbind
## Determine consensus peaks and create counts table

# Make a dba object from 
ATAC_dba <- dba(sampleSheet = ATAC_allnarrowpeaks)

# PCA plot
ATAC_dba_PCA <- dba.plotPCA(ATAC_dba, DBA_CONDITION)
plot(ATAC_dba)

# Write plots to .pngs
Cairo::Cairo(file ="20230508_NUDUL1_FOXO1_ATAC_diffbind_chrMrm_PCA_nonorm_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
ATAC_dba_PCA
dev.off()

Cairo::Cairo(file ="20230508_NUDUL1_FOXO1_ATAC_diffbind_chrMrm_corrheatmap_nonorm_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
plot(ATAC_dba)
dev.off()

# Determine consensus peaks, minOverlap = min # of peaksets a peak must be in to be included
ATAC_dba_consensus <- dba.peakset(ATAC_dba, consensus = DBA_CONDITION, minOverlap = 2)
ATAC_dba_consensus <- dba(ATAC_dba_consensus, mask = ATAC_dba_consensus$masks$Consensus, minOverlap = 1)

consensus_peaks <- dba.peakset(ATAC_dba_consensus, bRetrieve = TRUE)
consensus_peaks_df <- as.data.frame(consensus_peaks)

write.table(consensus_peaks_df, file = "20230508 NUDUL1_FOXO1_ATAC_consensuspeaks_chrMrm.txt", sep = "\t", row.names = FALSE, quote =FALSE)

# Create counts table from consensus_peaks
ATAC_dba$config$singleEnd <- FALSE
ATAC_dba_count <- dba.count(ATAC_dba, score = DBA_SCORE_READS, peaks = consensus_peaks, bRemoveDuplicates = FALSE, bUseSummarizeOverlaps = TRUE) 

# Create peakset with counts table
ATAC_peakset_counts <- dba.peakset(ATAC_dba_count, bRetrieve = TRUE)
ATAC_peakset_counts_df <- as.data.frame(ATAC_peakset_counts)

write.table(ATAC_peakset_counts_df, file = "20230508_NUDUL1_FOXO1_ATAC_peakset_counts_df.txt", sep = "\t", row.names = FALSE, quote =FALSE)


ATAC_dba_count_SE <- dba(ATAC_dba_count, bSummarizedExperiment = TRUE)

write.table(cbind(data.frame(chr = seqnames(ATAC_peakset_counts), start = start(ATAC_peakset_counts), end = end(ATAC_peakset_counts), strand = strand(ATAC_peakset_counts)), assay(ATAC_dba_count_SE)), file = "20230508_NUDUL1_FOXO1_ATAC_countstable_peakset_SE_chrMrm.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## Run DESEQ
### Normalize data
###Load counts table and process for DESeq
# Load counts table populated with peakset info (from line 89)
count <- read.table("20230508_NUDUL1_FOXO1_ATAC_countstable_peakset_SE_chrMrm.txt", sep = "\t", header = TRUE)
 
# Create DESeq Dataset 
countData <- count[ , c(5:14)]
colData <- data.frame(condition = factor(c("0hr", "0hr", "15min", "15min", "30min", "30min", "1hr", "1hr", "2hr", "2hr"), levels = c("0hr", "15min", "30min", "1hr", "2hr")), sample = c("0hrA", "0hrB", "15minA", "15minB", "30minA", "30minB", "1hrA", "1hrB", "2hrA", "2hrB"))

#### Load read counts
options(scipen = 99)
idxstats_01 <- readr::read_delim("idxstats_9330-HL-0001_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_02 <- readr::read_delim("idxstats_9330-HL-0002_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_03 <- readr::read_delim("idxstats_9330-HL-0003_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_04 <- readr::read_delim("idxstats_9330-HL-0004_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_05 <- readr::read_delim("idxstats_9330-HL-0005_S1_L005_samp81mill.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_06 <- readr::read_delim("idxstats_9330-HL-0006_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_07 <- readr::read_delim("idxstats_9330-HL-0007_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_08 <- readr::read_delim("idxstats_9330-HL-0008_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_09 <- readr::read_delim("idxstats_9330-HL-0009_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))
idxstats_10 <- readr::read_delim("idxstats_9330-HL-0010_S1_L005.hg19dm3.F4q10.sorted.txt", col_names = c("chr", "total", "actual", "n"))

idxstatsadj <- function(df, title) {
  df %>% 
    dplyr::mutate(sample = title) %>%
    dplyr::select(sample, chr, actual) %>%
    tidyr::pivot_wider(names_from = chr, values_from = actual)
} 

idxstats_01_adj <- idxstatsadj(idxstats_01, "9330-1")
idxstats_02_adj <- idxstatsadj(idxstats_02, "9330-2")
idxstats_03_adj <- idxstatsadj(idxstats_03, "9330-3")
idxstats_04_adj <- idxstatsadj(idxstats_04, "9330-4")
idxstats_05_adj <- idxstatsadj(idxstats_05, "9330-5")
idxstats_06_adj <- idxstatsadj(idxstats_06, "9330-6")
idxstats_07_adj <- idxstatsadj(idxstats_07, "9330-7")
idxstats_08_adj <- idxstatsadj(idxstats_08, "9330-8")
idxstats_09_adj <- idxstatsadj(idxstats_09, "9330-9")
idxstats_10_adj <- idxstatsadj(idxstats_10, "9330-10")

readcounts <- idxstats_01_adj %>%
  dplyr::add_row(idxstats_02_adj) %>%
  dplyr::add_row(idxstats_03_adj) %>%
  dplyr::add_row(idxstats_04_adj) %>%
  dplyr::add_row(idxstats_05_adj) %>%
  dplyr::add_row(idxstats_06_adj) %>%
  dplyr::add_row(idxstats_07_adj) %>%
  dplyr::add_row(idxstats_08_adj) %>%
  dplyr::add_row(idxstats_09_adj) %>%
  dplyr::add_row(idxstats_10_adj) %>%
  dplyr::mutate(hg19 = rowSums(.[2:26])) %>%
  dplyr::mutate(dm3 = rowSums(.[27:41])) %>%
  dplyr::select(sample, hg19, dm3)
readr::write_delim(readcounts, "20230508 NUDUL1_FOXO1_ATAC_readcounts_hg19_dm3.txt", delim = "\t")


#### Total Counts

# Create DESeq Dataset (cont'd)
dds_hg19chrMrm <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

#input the total counts from the .bam files using samtools idxstats (calculated in excel)
hg19chrMrm <- readcounts$hg19
sizeFactors(dds_hg19chrMrm) <- (hg19chrMrm/median(hg19chrMrm))

#print out normalization factors determined by DESeq2
sizeFactorsddsATAC_hg19chrMrm <- sizeFactors(dds_hg19chrMrm)
sizeFactorsddsATAC_hg19chrMrm

write.table(sizeFactorsddsATAC_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_ATAC_sizefactorsdds_hg19chrMrm.txt", sep = "\t", row.names = FALSE, quote =FALSE)

#dds <- DESeq(dds)
dds_hg19chrMrm <- estimateDispersions(dds_hg19chrMrm)
dds_hg19chrMrm <- nbinomWaldTest(dds_hg19chrMrm)
 
vsd <- varianceStabilizingTransformation(dds_hg19chrMrm)
PCA_hg19chrMrm <- plotPCA(vsd, intgroup = c("condition"))+ 
  geom_text(aes(label = vsd$sample), size=2) +
  ggtitle("ATAC ChIPseq PCA (hg19chrMrm norm)") +
  scale_color_manual(values = viridis(7)) +
  theme_classic()
PCA_hg19chrMrm


#### S2

# Create DESeq Dataset (cont'd)
dds_dm3<- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

#input the total counts from the .bam files using samtools idxstats (calculated in excel)
dm3<- readcounts$dm3
sizeFactors(dds_dm3) <- (dm3/median(dm3))

#print out normalization factors determined by DESeq2
sizeFactorsddsATAC_dm3<- sizeFactors(dds_hg19chrMrm)
sizeFactorsddsATAC_dm3

write.table(sizeFactorsddsATAC_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_ATAC_sizefactorsdds_dm3.txt", sep = "\t", row.names = FALSE, quote =FALSE)

#dds <- DESeq(dds)
dds_dm3<- estimateDispersions(dds_dm3)
dds_dm3<- nbinomWaldTest(dds_dm3)
 
vsd <- varianceStabilizingTransformation(dds_dm3)
PCA_dm3<- plotPCA(vsd, intgroup = c("condition"))+ 
  geom_text(aes(label = vsd$sample), size=2) +
  ggtitle("ATAC ChIPseq PCA (dm3norm)") +
  scale_color_manual(values = viridis(7)) +
  theme_classic()
PCA_dm3


#### DEseq Normalization
# Create DESeq Dataset (cont'd)
dds_deseq <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

#Estimate size factors with DESeq
dds_deseq <- estimateSizeFactors(dds_deseq)
 
#print out normalization factors determined by DESeq2
sizeFactorsddsATAC_deseq <- sizeFactors(dds_deseq)
sizeFactorsddsATAC_deseq

## Reciprocal size factors for deeptools 
print(1/sizeFactorsddsATAC_deseq) 

write.table(sizeFactorsddsATAC_deseq, file = "20230508_NUDUL1_FOXO1_ATAC_sizefactorsdds_deseq_hg19chrMrm.txt", sep = "\t", row.names = FALSE, quote =FALSE)

#dds <- DESeq(dds)
dds_deseq <- estimateDispersions(dds_deseq)
dds_deseq <- nbinomWaldTest(dds_deseq)
 
vsd_deseq <- varianceStabilizingTransformation(dds_deseq)
PCA_deseq <- plotPCA(vsd_deseq, intgroup = c("condition"))+ 
  geom_text(aes(label = vsd_deseq$sample), size=2) +
  ggtitle("ATAC ChIPseq PCA (deseq_hg19chrMrm norm)") +
  scale_color_manual(values = viridis(7)) +
  theme_classic()
PCA_deseq 


#### Print PCA plots to .png
Use grid.arrange to create one .png with all three plots

Cairo::Cairo(file ="20230508_NUDUL1_FOXO1_ATAC_diffbind_PCA_hg19chrMrm_deseq_5x7.5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 7.5, 
             pointsize = 10, 
             dpi = 300)
grid.arrange(PCA_hg19chrMrm, PCA_dm3, PCA_deseq)
dev.off()


### Calculate differential changes uing DESeq normalization
normdata <- counts(dds_deseq, normalized = TRUE)
colnames(normdata) <- paste("norm",c("0hrA", "0hrB", "15minA", "15minB", "30minA", "30minB", "1hrA", "1hrB", "2hrA", "2hrB"), sep = "_")

result_dTAG47_15minvs0hr <- results(dds_deseq, contrast = c("condition","15min","0hr"))
ATAC_15v0_count_norm_results_deseq <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_15minvs0hr))
write.table(ATAC_15v0_count_norm_results_deseq, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_15minv0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_30minvs0hr <- results(dds_deseq, contrast = c("condition","30min","0hr"))
ATAC_30v0_count_norm_results_deseq <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_30minvs0hr))
write.table(ATAC_30v0_count_norm_results_deseq, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_30v0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_1hrvs0hr <- results(dds_deseq, contrast = c("condition","1hr","0hr"))
ATAC_1v0_count_norm_results_deseq <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_1hrvs0hr))
write.table(ATAC_1v0_count_norm_results_deseq, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_1hrv0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_2hrvs0hr <- results(dds_deseq, contrast = c("condition","2hr","0hr"))
ATAC_2v0_count_norm_results_deseq <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_2hrvs0hr))
write.table(ATAC_2v0_count_norm_results_deseq, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_deseqnorm_hg19chrMrm_2hrv0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")


### Make MA plots
Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_deseq_hg19chrMrm_0v15min_minus7.5_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_15minvs0hr, ylim = c(-3, 3), main = "ATAC 15min (+dTAG) deseq", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_deseq_hg19chrMrm_0v30min_minus7.5_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_30minvs0hr, ylim = c(-3, 3), main = "ATAC 30min (+dTAG) deseq", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_deseq_hg19chrMrm_0v1hr_minus7.5_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_1hrvs0hr, ylim = c(-3, 3), main = "ATAC 1 hr (+dTAG) deseq", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_deseq_hg19chrMrm_0v2hr_minus7.5_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_2hrvs0hr, ylim = c(-3, 3), main = "ATAC 2 hr (+dTAG) deseq", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

### Calculate differential changes uing total counts normalization
normdata_hg19chrMrm <- counts(dds_hg19chrMrm, normalized = TRUE)
colnames(normdata_hg19chrMrm) <- paste("norm",c("0hrA", "0hrB", "15minA", "15minB", "30minA", "30minB", "1hrA", "1hrB", "2hrA", "2hrB"), sep = "_")
 
result_dTAG47_15minvs0hr_hg19chrMrm <- results(dds_hg19chrMrm, contrast = c("condition","15min","0hr"))
ATAC_15v0_count_norm_results_hg19chrMrm <- as.data.frame(cbind(count[,c(1:3)], countData, normdata_hg19chrMrm, result_dTAG47_15minvs0hr_hg19chrMrm))
write.table(ATAC_15v0_count_norm_results_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_hg19chrMrm_15v0hr_minus7.5.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_30minvs0hr_hg19chrMrm <- results(dds_hg19chrMrm, contrast = c("condition","30min","0hr"))
ATAC_30v0_count_norm_results_hg19chrMrm <- as.data.frame(cbind(count[,c(1:3)], countData, normdata_hg19chrMrm, result_dTAG47_30minvs0hr_hg19chrMrm))
write.table(ATAC_30v0_count_norm_results_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_hg19chrMrm_30v0hr_minus7.5.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_1hrvs0hr_hg19chrMrm <- results(dds_hg19chrMrm, contrast = c("condition","1hr","0hr"))
ATAC_1v0_count_norm_results_hg19chrMrm <- as.data.frame(cbind(count[,c(1:3)], countData, normdata_hg19chrMrm, result_dTAG47_1hrvs0hr_hg19chrMrm))
write.table(ATAC_1v0_count_norm_results_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_hg19chrMrm_1v0hr_minus7.5.txt", quote = FALSE, row.names = FALSE, sep = "\t")

result_dTAG47_2hrvs0hr_hg19chrMrm <- results(dds_hg19chrMrm, contrast = c("condition","2hr","0hr"))
ATAC_2v0_count_norm_results_hg19chrMrm <- as.data.frame(cbind(count[,c(1:3)], countData, normdata_hg19chrMrm, result_dTAG47_2hrvs0hr_hg19chrMrm))
write.table(ATAC_2v0_count_norm_results_hg19chrMrm, file = "20230508 NUDUL1_FOXO1_FKBP_ATAC_hg19chrMrm_2v0hr_minus7.5.txt", quote = FALSE, row.names = FALSE, sep = "\t")

### Make MA plots, total counts
Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_hg19norm_hg19chrMrm_0v15_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_15minvs0hr_hg19chrMrm, ylim = c(-3, 3), main = "ATAC 15min (+dTAG) hg19chrMrm hg19 norm", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_hg19norm_hg19chrMrm_0v30_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_30minvs0hr_hg19chrMrm, ylim = c(-3, 3), main = "ATAC 30min (+dTAG) hg19chrMrm hg19 norm", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_hg19norm_hg19chrMrm_0v1_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_1hrvs0hr_hg19chrMrm, ylim = c(-3, 3), main = "ATAC 1 hr (+dTAG) hg19chrMrm hg19 norm", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

Cairo::Cairo(file ="20230508 NUDUL1_FOXO1_ATAC_hg19norm_hg19chrMrm_0v2_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_2hrvs0hr_hg19chrMrm, ylim = c(-3, 3), main = "ATAC 2 hr (+dTAG) hg19chrMrm hg19 norm", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

## End
```
## Creation of bigwig files for visualization.
Bigwig files were created using DeepTools in conda. Files were normalized using the reciprocal of the DESeq sizeFactor.
```
bamCoverage --bam /Volumes/Hillary_2p/9330-NUDUL1_ATAC/9330-HL-0003_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.bam -o /Volumes/Hillary_2p/9330-NUDUL1_ATAC/9330-HL-0003_S1_L005.hg19dm3.F4q10.sorted.chrMrm.sorted.sfdeseqrecip.bs10.bw --binSize 10 --scaleFactor 1.3087450 --effectiveGenomeSize 2827437033 --blackListFileName /Volumes/Hillary_2p/genomes/hg19_blacklist/hg19.blacklistpeaks.bed --extendReads\
```
