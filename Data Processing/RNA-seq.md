# RNA-seq Data Processing
The code used to process RNA-seq data is presented below.
## Adaptor Trimming
Adaptor trimming was performed using Trimmomatic in bash. Please note that file names may need to be changed based on naming conventions at the time of sequencing.
```
for file in $(ls *_R1_001.fastq.gz); #Lis/Users/kristystengel/Downloads/Trimmomatic_PE_KSts all R1.fastq files in directory
do
base=${file//_R1_001.fastq.gz} #Extracts base filename, ie XXXX_HL-X_S1_L005
if [ ! -f "${base}"-forward_paired.fastq.gz ]; #Executes following code if there is no corresponding index file
then
echo "${base}" #Prints base
java -classpath /Users/jing/software/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 "${base}"_R1_001.fastq.gz "${base}"_R2_001.fastq.gz "${base}"-forward_paired.fastq.gz "${base}"-forward_unpaired.fastq.gz "${base}"-reverse_paired.fastq.gz "${base}"-reverse_unpaired.fastq.gz ILLUMINACLIP:../Hillary/TruSeq_CD_adapter.txt:2:30:10:1:TRUE LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 #Trims files
echo "Trimming of "${base}"_R1_001.fastq.gz "${base}"_R2_001.fastq.gz complete" #Prints string to console when trimming is complete
else
echo "Skipping "${base}" as "${base}"-forward_paired.fastq.gz already exists" #Prints string to console if corresponding BAM file already exists
fi
done
```
## Alignment
Trimmed data was aligned to human (hg19)(dm3) genome using Tophat in bash.
```
for file in $(ls *-forward_paired.fastq.gz); #Lists all R1.fastq files in directory
do
if [ ! -f ${file//-forward_paired.fastq.gz}.hg19dm3.sam ]; #Executes following code if SAM file doesn't exist
then
base=${file//-forward_paired.fastq.gz} #Extracts base filename, ie XXXX_HL-X_S1_L005
echo "${base}" #Prints base
tophat -p 8 -G /Users/yuezhao/DATA1/seq-scripts/UCSC-annotations/Homo_sapiens-hg19/UCSC/hg19/Annotation/Genes/UCSC-hg19-genes.gtf -o "${base}"-tophat /Users/yuezhao/DATA1/seq-scripts/UCSC-annotations/Homo_sapiens-hg19/UCSC/hg19/Sequence/Bowtie2Index/genome "${base}"-forward_paired.fastq.gz "${base}"-reverse_paired.fastq.gz
echo "Tophat alignment of -1 "${base}"-forward_paired.fastq.gz -2 "${base}"-reverse_paired.fastq.gz to hg19dm3 complete" #Prints string to console when trimming is complete
else
echo "Skipping alignment of "${file//-forward_paired.fastq.gz}" as "${file//-forward_paired.fastq.gz}".hg19.sam already exists" #Prints string to console if SAM file already exists
fi
done
```

## Transcript quantification and differential expression 
Quantification of transcripts and differential expression was performed using CuffDiff in bash. 
```
NUDUL1_par=("9718-HL-0001_S1_L005-tophat/accepted_hits.bam,9718-HL-0002_S1_L005-tophat/accepted_hits.bam")
NUDUL1_0h=("9718-HL-0003_S1_L005-tophat/accepted_hits.bam,9718-HL-0004_S1_L005-tophat/accepted_hits.bam")
cuffdiff -o NUDUL1_parvfkbp/ -b /Users/hiebertlab/DATA1/seq-scripts/UCSC-annotations/Homo_sapiens-hg19/UCSC/hg19/Sequence/Bowtie2Index/genome.fa -p 8 -L par,fkbp -u /Users/hiebertlab/DATA1/seq-scripts/UCSC-annotations/Homo_sapiens-hg19/UCSC/hg19/Annotation/Genes/UCSC-hg19-genes.gtf ${NUDUL1_par} ${NUDUL1_0h}
```
