# PRO-seq Data Processing
The code used to process PRO-seq data is presented below.
## Adaptor Trimming
Adaptor trimming was performed using Trimmomatic in bash. Please note that file names may need to be changed based on naming conventions at the time of sequencing.
```
for file in $(ls *.fastq.gz); #Lists all .fastq files in directory
do
java -classpath /Users/jing/software/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 4 $file ${file%.*}.trimmed.txt CROP:75 ILLUMINACLIP:../adaptor-8.txt:2:30:7 TRAILING:15 MINLEN:15 #Trims files
echo "Trimming of $file to ${file%.*}.trimmed.txt complete" #Prints string to console when trimming is complete
done
```
## Create reverse complement
Trimmed data was converted to the reverse complement using the Fastx toolkit in bash.
```
for file in $(ls *.txt); #Creates a list of .txt files in working directory
do
if [ ! -f ${file%.*}.revcomp.txt ]; #Executes following code if .txt does not have corresponding revcomp file
then
echo "Creating reverse complement for $file" #Prints string to console
fastx_reverse_complement -Q33 -i $file -o ${file%.*}.revcomp.txt
else
echo "Skipping $file as ${file%.*}.revcomp.txt already exists" #Prints string to console
fi
done
```

## Alignment
Reverse complemented data was aligned to a concatenated human (hg19)/ Drosophila (dm3) genome using Bowtie2 in bash.
```
for file in $(ls *.revcomp.txt); #Creates a .revcomp.txt files in working directory
do
if [ ! -f $file_hg19dm3.sam ]; #Executes following code if .txt does not have corresponding sam file
then
echo "Aligning $file to hg19_dm3" #Prints string to console
bowtie2 -p 4 --mm -x /Volumes/Hillary_2p/genomes/hg19_dm3/hg19+dm3 -U $file -S ${file%.*}_hg19dm3.sam
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
## Nascent transcript and eRAN quantification and differential expression 
Quantification of  gene body and intergenic polymerases and differential expression was performed using NRSA in bash. 
```
## Gene body
nohup perl /Users/jing/software/NRSA-v2/bin/pause_PROseq.pl -o LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &nohup perl /Users/jing/software/NRSA-v2/bin/pause_PROseq.pl -o LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &nohup perl /Users/jing/software/NRSA-v2/bin/pause_PROseq.pl -o LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &nohup perl /Users/jing/software/NRSA-v2/bin/pause_PROseq.pl -o LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &nohup perl /Users/jing/software/NRSA-v2/bin/pause_PROseq.pl -o LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &

## eRNA
nohup perl /Users/jing/software/NRSA-v2/bin/eRNA.pl -w LY1_FOXO1_0hrv2_DEseq/ -in1 7771-HL-1_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam -in2 7771-HL-7_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam 7771-HL-8_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_0v2.out &
```

## Creation of bedgraph files for visualization.
Bedgraph files were created using Homer in bash. 
```
## Make tag directories
nohup makeTagDirectory LY1_FOXO1_PS_0hrB/ 7771-HL-2_S1_L005_R1_001.fastq.trimmed.revcomp_hg19dm3.F4q10.sorted.dm3rm.sorted.bam > nohup_02.out &

## Make stranded bedgraphs
for file in $(ls -d */); #Creates a list of tag directories
do
if [ ! -f ${file%/*}.plus.bedgraph.gz ]; #Executes following code if there is no corresponding bedgraph file
then
echo "Making plus strand UCSC file from $file" #Prints string to console
makeUCSCfile $file -strand + -fragLength 66 -o ${file%/*}.plus.bedgraph
else
echo "Skipping $file as ${file%/*}.plus.bedgraph already exists" #Prints string to console if corresponding bedgraph file already exists
fi
if [ ! -f ${file%/*}.minus.bedgraph.gz ]; #Executes following code if there is no corresponding bedgraph file
then
echo "Making minus strand UCSC file from $file" #Prints string to console
makeUCSCfile $file -strand - -neg -fragLength 66 -o ${file%/*}.minus.bedgraph
else
echo "Skipping $file as ${file%/*}.minus.bedgraph already exists" #Prints string to console if corresponding bedgraph file already exists
fi
done
```
