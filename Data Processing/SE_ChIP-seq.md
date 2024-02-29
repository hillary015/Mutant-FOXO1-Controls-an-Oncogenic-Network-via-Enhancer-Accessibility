# ChIP-seq Data Processing
The code used to process single-end ChIP-seq (H3K4me1) data is presented below.
## Adaptor Trimming
Adaptor trimming was performed using Trimmomatic in bash. Please note that file names may need to be changed based on naming conventions at the time of sequencing.
```
for file in $(ls *.fastq.gz); #Lists all .fastq files in directory
do
java -classpath /Users/hiebertlab/software/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 8 $file ${file%.*}.trimmed.txt ILLUMINACLIP:../TruSeq_CD_adapter.txt:2:30:10:1 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 #Trims files
echo "Trimming of $file to ${file%.*}.trimmed.txt complete" #Prints string to console when trimming is complete
done
```
## Alignment
Trimmed data was aligned to a concatenated human (hg19)/ Drosophila (dm3) genome using Bowtie2 in bash.
```
for file in $(ls *.trimmed.txt); #Creates a .revcomp.txt files in working directory
do
if [ ! -f $file_hg19dm3.sam ]; #Executes following code if .txt does not have corresponding sam file
then
echo "Aligning $file to hg19_dm3" #Prints string to console
bowtie2 -p 8 --mm -x /Volumes/Hillary_3/genomes/hg19_dm3/hg19+dm3 -U $file -S ${file%.*}_hg19dm3.sam
else
echo "Skipping $file as ${file%.*}_hg19dm3.sam already exists" #Prints string to console
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
## Creation of bigwig files for visualization.
Bigwig files were created using DeepTools in conda.
```
bamCoverage --bam LY1_H3K4me1_NT.fastq.trimmed_hg19dm3.F4q10.sorted.chrMrm.sorted.bam -o LY1_H3K4me1_NT.fastq.trimmed_hg19dm3.F4q10.sorted.chrMrm.sorted.sfdeseqrecip.bs10.bw --binSize 10 --effectiveGenomeSize 2827437033 --blackListFileName /Volumes/Hillary_X6/genomes/hg19_blacklist/hg19.blacklistpeaks.bed --extendReads 185
