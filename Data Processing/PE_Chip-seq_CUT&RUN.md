# ChIP-seq Data Processing
The code used to process single-end ChIP-seq (H3K4me1) data is presented below.
## Adaptor Trimming
Adaptor trimming was performed using Trimmomatic in bash. Please note that file names may need to be changed based on naming conventions at the time of sequencing.
```
nohup java -classpath /Users/yuezhao/Software/Trimmomatic-0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 5400-HL-1-GAGAGTAC-TGGAAGCA_S1_L001_R1_001.fastq.gz 5400-HL-1-GAGAGTAC-TGGAAGCA_S1_L001_R2_001.fastq.gz 5400-HL-1-forward_paired.fastq.gz 5400-HL-1-forward_unpaired.fastq.gz 5400-HL-1-reverse_paired.fastq.gz 5400-HL-1-reverse_unpaired.fastq.gz ILLUMINACLIP:../../TruSeq_CD_adapter.txt:2:30:10:1:TRUE LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 > nohup-1.out &
```
## Alignment
Trimmed data was aligned to a concatenated human (hg19)/ Drosophila (dm3) genome (ChIP) or human (hg19)/ E. coli (sc) (CUT&RUN) genome using Bowtie2 in bash.
```
nohup bowtie2 -p 8 -x /Volumes/Clare/Clare/hg19+dm3 -1 5400-HL-8-forward_paired.fastq.gz -2 5400-HL-8-reverse_paired.fastq.gz --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -S 5400-HL-8_LY1_H3K27ac_24hrB.hg19dm.sam >> nohup-08.out &
```
## SAM/BAM conversion and read filtering
Aligned SAM files were converted to BAM files and low-quality reads were removed using SAMtools in bash.
```
for file in $(ls *.sam); #Creates a list of SAM files in working directory
do
if [ ! -f ${file%.*}.bam ]; #Executes following code if SAM does not have corresponding BAM
then
echo "Converting $file to ${file%.*}.bam" #Prints string to console
samtools view -S -b -@16 $file > ${file%.*}.bam # Converts SAM file to BAM file
else
echo "Skipping conversion of $file as ${file%.*}.bam already exists" #Prints string to console if corresponding BAM file already exists
fi

if [ ! -f ${file%.*}.sorted.bam ]; #Executes the following code if BAM file hasn't been sorted
then
echo "Sorting ${file%.*}.bam and creating new file ${file%.*}.sorted.bam" #Prints string to console
samtools sort -@16 -T temp ${file%.*}.bam -o ${file%.*}.sorted.bam  #Resorts BAM file
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
## Count and remove spikein
```
# Count reads aligned to each chromosome, extract hg19 reads, resort and index new bam files
for file in $(ls *.F4q10.sorted.bam); #Creates a list of sorted BAM files in working directory
do
if [ ! -f idxstats_${file%.*}.txt ]; #Executes following code if sorted.BAM does not have corresponding txt file
then
echo "Reporting alignment statistics for $file and printing to idxstats_${file%.*}.txt" #Prints string to console
samtools idxstats $file > idxstats_${file%.*}.txt #Reports alignment statistics for sorted.bam file and prints to .txt file
else
echo "Skipping $file as ${file%.*}.sorted.txt already exists" #Prints string to console if corresponding BAM file already exists
fi

if [ ! -f ${file%.*}.chrMrm.bam ]; #Executes following code if mito reads have not yet been removed from sorted.bam file
then
echo "Printing hg19 reads from $file to ${file%.*}.chrMrm.bam" #Prints string to console
samtools view -bh $file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file%.*}.chrMrm.bam
else
echo "Skipping $file as ${file%.*}.chrMrm.bam already exists" #Prints string to console if corresponding BAM file already exists
fi

if [ ! -f ${file%.*}.chrMrm.sorted.bam ]; #Executes following code if BAM file hasn't been resorted
then
echo "Sorting ${file%.*}.bam and creating new file ${file%.*}.sorted.bam" #Prints string to console
samtools sort -@ 16 -T ${file%.*}.temp -o ${file%.*}.chrMrm.sorted.bam ${file%.*}.chrMrm.bam #Reorts BAM file
else
echo "Skipping sorting of ${file%.*}.bam as ${file%.*}.sorted.bam already exists" #Prints string to console if resorted BAM file already exists
fi

if [ ! -f ${file%.*}.sorted.bam.bai ]; #Executes following code if BAM file hasn't been indexed
then
echo "Indexing ${file%.*}.sorted.bam and creating new file ${file%.*}.sorted.bam.bai" #Prints string to console
samtools index ${file%.*}.chrMrm.sorted.bam #Indexes BAM file
else
echo "Skipping indexing of ${file%.*}.sorted.bam as ${file%.*}.sorted.bam.bai already exists" #Prints string to console if indexed BAM file already exists
fi
done
```
## Call peaks with MACS2
Peaks were called in using MACS2 in conda.
```
macs2 callpeak -B -t DHL4_BCL6.fastq.trimmed_hg19dm3.F4q10.sorted.chrMrm.sorted.bam -c DHL4_input_NT.fastq.trimmed_hg19dm3.F4q10.sorted.chrMrm.sorted.bam -f BAM -g hs -q 0.01 --SPMR --call-summits --outdir macs2 -n DHL4_BCL6_NT.fastq.trimmed_hg19dm3.F4q10.sorted.chrMrm.sorted.callpeaks.q0.01.SPMR.callsummits
```
## Quantification of reads and differential expression analysis
A consensus peakset was defined using DiffBind in R. Reads within the peaks defined in the consensus peakset were quantified with DiffBind in R. Counts were normalized and differential expression analysis was performed with DESeq2 in R.
```
## Load libraries
library("DiffBind")
packageVersion("DiffBind")
library("DESeq2")
packageVersion("DESeq2")
library(tidyverse)
library(gridExtra)
library(viridi

## Load Data
Data loaded from external hard drive (.bai files corresponding to .bam files loaded need to be located in the same directory), loaded into samplesheet format
H3K27ac_allnarrowpeaks <- data.frame(SampleID = c("0hrA", "0hrB", "6hrA", "6hrB"),
                                  Condition = c("0hr", "0hr", "6hr", "6hr"),
                                  Replicate = c("A", "B", "A", "B"),
                                  bamReads = c("D:/5400-LY1_H3K27ac/5400-HL-1_LY1_H3K27ac_0hrA.hg19dm.F4q10.sorted.dmrm.sorted.bam", "D:/5400-LY1_H3K27ac/5400-HL-2_LY1_H3K27ac_0hrB.hg19dm.F4q10.sorted.dmrm.sorted.bam", "D:/5400-LY1_H3K27ac/5400-HL-5_LY1_H3K27ac_6hrA.hg19dm.F4q10.sorted.dmrm.sorted.bam", "D:/5400-LY1_H3K27ac/5400-HL-6_LY1_H3K27ac_6hrB.hg19dm.F4q10.sorted.dmrm.sorted.bam"),
                                  Peaks = c("D:/5400-LY1_H3K27ac/5400-HL-1_LY1_H3K27ac_0hrA.hg19dm.F4q10.sorted.dmrm.sorted.callpeaks.q0.01.SPMR.callsummits_peaks.narrowPeak", "D:/5400-LY1_H3K27ac/5400-HL-2_LY1_H3K27ac_0hrB.hg19dm.F4q10.sorted.dmrm.sorted.callpeaks.q0.01.SPMR.callsummits_peaks.narrowPeak",  "D:/5400-LY1_H3K27ac/5400-HL-5_LY1_H3K27ac_6hrA.hg19dm.F4q10.sorted.dmrm.sorted.callpeaks.q0.01.SPMR.callsummits_peaks.narrowPeak", "D:/5400-LY1_H3K27ac/5400-HL-6_LY1_H3K27ac_6hrB.hg19dm.F4q10.sorted.dmrm.sorted.callpeaks.q0.01.SPMR.callsummits_peaks.narrowPeak"), PeakCaller = "narrow")

## Diffbind
Determine consensus peaks and create counts table
# Make a dba object from 
H3K27ac_dba <- dba(sampleSheet = H3K27ac_allnarrowpeaks)

# PCA plot
H3K27ac_dba_PCA <- dba.plotPCA(H3K27ac_dba, DBA_CONDITION)
plot(H3K27ac_dba)

# Write plots to .pngs
Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_diffbind_PCA_nonorm_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
H3K27ac_dba_PCA
dev.off()

Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_diffbind_corrheatmap_nonorm_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 12, 
             dpi = 300)
plot(H3K27ac_dba)
dev.off()

# Determine consensus peaks, minOverlap = min # of peaksets a peak must be in to be included
H3K27ac_dba_consensus <- dba.peakset(H3K27ac_dba, consensus = DBA_CONDITION, minOverlap = 1)
H3K27ac_dba_consensus <- dba(H3K27ac_dba_consensus, mask = H3K27ac_dba_consensus$masks$Consensus, minOverlap = 1)

consensus_peaks <- dba.peakset(H3K27ac_dba_consensus, bRetrieve = TRUE)
consensus_peaks_df <- as.data.frame(consensus_peaks)

write.table(consensus_peaks_df, file = "20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_consensuspeaks_mo1.txt", sep = "\t", row.names = FALSE, quote =FALSE)

# Create counts table from consensus_peaks
H3K27ac_dba$config$singleEnd <- FALSE
H3K27ac_dba_count <- dba.count(H3K27ac_dba, score = DBA_SCORE_READS, peaks = consensus_peaks, bRemoveDuplicates = FALSE, bUseSummarizeOverlaps = TRUE) 

# Create peakset with counts table
H3K27ac_peakset_counts <- dba.peakset(H3K27ac_dba_count, bRetrieve = TRUE)
H3K27ac_peakset_counts_df <- as.data.frame(H3K27ac_peakset_counts)

write.table(H3K27ac_peakset_counts_df, file = "20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_peakset_counts_df_mo1.txt", sep = "\t", row.names = FALSE, quote =FALSE)


H3K27_dba_count_SE <- dba(H3K27ac_dba_count, bSummarizedExperiment = TRUE)

write.table(cbind(data.frame(chr = seqnames(H3K27ac_peakset_counts), start = start(H3K27ac_peakset_counts), end = end(H3K27ac_peakset_counts), strand = strand(H3K27ac_peakset_counts)), assay(H3K27_dba_count_SE)), file = "20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_coutstable_peakset_mo1_SE.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## Run DESEQ
### Normalize data
Load counts table and process for DESeq
```{r}
# Load counts table populated with peakset info (from line 89)
count <- read.table("20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_coutstable_peakset_mo1_SE.txt", sep = "\t", header = TRUE)
 
# Create DESeq Dataset 
countData <- count[ ,c(5:8)]
colData <- data.frame(condition = factor(c("0hr", "0hr", "6hr", "6hr"), levels = c("0hr", "6hr")), sample = c("0hrA", "0hrB", "6hrA", "6hrB"))

#### DEseq Normalization

# Create DESeq Dataset (cont'd)
dds_deseq <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

#Estimate size factors with DESeq
dds_deseq <- estimateSizeFactors(dds_deseq)
 
#print out normalization factors determined by DESeq2
sizeFactorsddsH3K27ac_deseq <- sizeFactors(dds_deseq)
sizeFactorsddsH3K27ac_deseq

write.table(sizeFactorsddsH3K27ac_deseq, file = "20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_sizefactorsdds_deseq.txt", sep = "\t", row.names = FALSE, quote =FALSE)

#dds <- DESeq(dds)
dds_deseq <- estimateDispersions(dds_deseq)
dds_deseq <- nbinomWaldTest(dds_deseq)
 
vsd_deseq <- varianceStabilizingTransformation(dds_deseq)
PCA_deseq <- plotPCA(vsd_deseq, intgroup = c("condition"))+ 
  geom_text(aes(label = vsd_deseq$sample), size=2) +
  ggtitle("H3K27ac ChIPseq PCA (deseq norm)") +
  scale_color_manual(values = viridis(5)) +
  theme_classic()
PCA_deseq 


#### Print PCA plots to .png
Use grid.arrange to create one .png with all three plots

Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_diffbind_PCA_tc_sc_deseq_5x7.5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 7.5, 
             pointsize = 10, 
             dpi = 300)
grid.arrange(PCA_tc, PCA_s2, PCA_deseq)
dev.off()

### Calculate differential changes uing DESeq normalization

normdata <- counts(dds_deseq, normalized = TRUE)
colnames(normdata) <- paste("norm",c("0hrA","0hrB","6hrA", "6hrB"), sep = "_")
 
result_dTAG47_6hrvs0hr <- results(dds_deseq, contrast = c("condition","6hr","0hr"))
H3K27ac_6v0_count_norm_results_deseq <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_6hrvs0hr))
write.table(H3K27ac_6v0_count_norm_results_deseq, file = "20240116LY1_FOXO1_FKBP_H3K27ac_deseqnorm_6v0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")


### Make MA plots

Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_H3K27ac_deseq_0v6_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_6hrvs0hr, ylim = c(-2.5, 2.5), main = "H3K27ac 6h (+dTAG) deseq", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()


### Calculate differential changes using spikein normalization

normdata_dm3 <- counts(dds_dm3, normalized = TRUE)
colnames(normdata_dm3) <- paste("norm",c("0hrA","0hrB", "6hrA", "6hrB"), sep = "_")
 
result_dTAG47_6hrvs0hr_dm3 <- results(dds_dm3, contrast = c("condition","6hr","0hr"))
H3K27ac_6v0_count_norm_results_dm3 <- as.data.frame(cbind(count[,c(1:3)], countData, normdata, result_dTAG47_6hrvs0hr_dm3))
write.table(H3K27ac_6v0_count_norm_results_dm3, file = "20240116_LY1_FOXO1_H3K27ac_0v6_dm3norm_6v0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")

### Make MA plots, spikein

Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_dm3_0v6_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_6hrvs0hr_dm3, ylim = c(-2.5, 2.5), main = "H3K27ac 6h (+dTAG) dm3", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()

### Calculate differential changes uing total counts normalization

normdata_tc <- counts(dds_tc, normalized = TRUE)
colnames(normdata_tc) <- paste("norm",c("0hrA","0hrB","6hrA", "6hrB"), sep = "_")
 
result_dTAG47_6hrvs0hr_tc <- results(dds_tc, contrast = c("condition","6hr","0hr"))
H3K27ac_6v0_count_norm_results_tc <- as.data.frame(cbind(count[,c(1:3)], countData, normdata_tc, result_dTAG47_6hrvs0hr_tc))
write.table(H3K27ac_6v0_count_norm_results_tc, file = "20240116_LY1_FOXO1_H3K27ac_0v6_tcnorm_6v0hr.txt", quote = FALSE, row.names = FALSE, sep = "\t")


### Make MA plots, total counts

Cairo::Cairo(file ="20240116_LY1_FOXO1_H3K27ac_0v6_tc_0v6_MA_5x5.png", 
             bg = "white",
             type = "png",
             units = "in", 
             width = 5, 
             height = 5, 
             pointsize = 10, 
             dpi = 300)
plotMA(result_dTAG47_6hrvs0hr_tc, ylim = c(-2.5, 2.5), main = "H3K27ac 6h (+dTAG) tc", cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75)
dev.off()
```
## Creation of bigwig files for visualization.
Bigwig files were created using DeepTools in conda. Files were normalized using the reciprocal of the DESeq sizeFactor.
```
bamCoverage --bam /Volumes/Clare/Hillary/5400-LY1-FOXO1_FKBP_H3K27ac_timecourse_ChIPseq_30millionreads/5400-HL-1_LY1_H3K27ac_0hrA.hg19dm.F4q10.sorted.dmrm.sorted.bam -o 5400-HL-1_LY1_H3K27ac_0hrA.hg19dm3.F4q10.sorted.dmrm.sorted.sfdeseqrecip.bs10.bw --binSize 10 --scaleFactor 2.443060207 --effectiveGenomeSize 2827437033 --blackListFileName /Volumes/Hillary_X6/genomes/hg19_blacklist/hg19.blacklistpeaks.bed --extendReads
```
