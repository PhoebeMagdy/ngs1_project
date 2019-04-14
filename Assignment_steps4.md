# Files Manipulation 

## 1. Download the files from SRA database
##```wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra```

## 2. Download to the SRAtoolkit
##```wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz```

## 3. unzip the SRAtoolkit
##```tar -xzf sratoolkit.2.9.6-ubuntu64.tar.gz```

## 4. Add SRAtoolkit to the PATH:
```PATH=$PATH: /home/ngs-01/Downloads/sratoolkit.2.9.6-ubuntu64/bin```

## 5. Convert *.SRA file to 2 separate FASTQ files (forward & reverse reads) using fastq-dump and downloading the first 5000000 reads only:
```fastq-dump --split-files -X 5000000 SRR8797509.sra```

## 6-counting number of sequences in the fastq file 
##```grep '@' SRR8797509.fastq | wc -l```

## 7-Use the seqkit tool to split the main fastq file into separate files
```seqkit split2  SRR8797509_1.fastq -s 1000000```
```seqkit split2  SRR8797509_2.fastq -s 1000000```

## 8- shuffleing the main fastq files (Forward and reverse reads) each containing the 5000000 reads
```seqkit shuffle SRR8797509_1.fastq > SRR8797509_1_shuffled.fastq```
```seqkit shuffle SRR8797509_2.fastq > SRR8797509_2_shuffled.fastq```

## 9-split the shuffled files 
```seqkit split2  SRR8797509_1_shuffled.fastq -s 1000000```
```seqkit split2  SRR8797509_2_shuffled.fastq -s 1000000```

# Quality Control checking

## checking the first sample before and after shuffling for quality control using FASTQC
## install the required tools for performing the quality control assessment
##```conda install -c bioconda fastqc```
##```conda install -c bioconda multiqc```

## FastQC assessment for the first sample before shuffling
```mkdir ~/workdir/Assignment/FASTQC_results && cd ~/workdir/Assignment/FASTQC_results```
```for f in ~/workdir/Assignment/SRR8797509.fastq.split/SRR8797509_*.part_001.fastq;do```
```fastqc -t 1 -f fastq -noextract -o . $f;done```

## merging the fastQC reports of the two reads of sample 1 before shuffling
```multiqc -z -o . .```

## FastQC assessment for the first sample after shuffling
mkdir ~/workdir/Assignment/FASTQC_results_shuffled && cd ~/workdir/Assignment/FASTQC_results_shuffled
for f in ~/workdir/Assignment/shuff_SRR8797509.fastq.split/shuff_SRR8797509_*.part_001.part_001.fastq;do 
fastqc -t 1 -f fastq -noextract -o . $f;done

## merging the fastQC reports of the two reads of sample 1 after shuffling
multiqc -z -o . .

# Data Trimming

## install the required tools for trimming data
## ```conda install -c bioconda trimmomatic```

## Mild trimming for the 5 unshuffled samples (5-forward && 5-reverse reads)


```javascript
mkdir ~/workdir/Assignment/Mild_trimmed && cd ~/workdir/Assignment/Mild_trimmed
for SAMPLE in 1 2 3 4 5;
    do
        R1=$HOME/workdir/Assignment/SRR8797509.fastq.split/SRR8797509_1.part_00${SAMPLE}.fastq
        R2=$HOME/workdir/Assignment/SRR8797509.fastq.split/SRR8797509_2.part_00${SAMPLE}.fastq
        newf1=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_1.part_00${SAMPLE}.pe.trim.fastq
        newf2=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_2.part_00${SAMPLE}.pe.trim.fastq
        newf1U=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_1.part_00${SAMPLE}.se.trim.fastq
        newf2U=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_2.part_00${SAMPLE}.se.trim.fastq

        adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

        trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $R1 $R2 $newf1 $newf1U $newf2 $newf2U \
        ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36 
    done
```


## Aggressive trimming for the 5 shuffled samples (5-forward && 5-reverse reads)
mkdir ~/workdir/Assignment/Aggressive_trimmed && cd ~/workdir/Assignment/Aggressive_trimmed
for SAMPLE in 1 2 3 4 5;
    do
        R1=$HOME/workdir/Assignment/shuff_SRR8797509.fastq.split/shuff_SRR8797509_1.part_001.part_00${SAMPLE}.fastq
        R2=$HOME/workdir/Assignment/shuff_SRR8797509.fastq.split/shuff_SRR8797509_2.part_001.part_00${SAMPLE}.fastq
        newf1=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_1.part_001.part_00${SAMPLE}.pe.trim.fastq
        newf2=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_2.part_001.part_00${SAMPLE}.pe.trim.fastq
        newf1U=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_1.part_001.part_00${SAMPLE}.se.trim.fastq
        newf2U=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_2.part_001.part_00${SAMPLE}.se.trim.fastq

        adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

        trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $R1 $R2 $newf1 $newf1U $newf2 $newf2U \
        ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:3:30 MINLEN:35 
    done

# Data Alignment

## install the reqired tool to make alignment with BWA
##```conda install -c bioconda bwa ```

## indexing 
```mkdir ~/workdir/Assignment/bwa_align/bwaIndex && cd ~/workdir/Assignment/bwa_align/bwaIndex```
```ln -s ~/workdir/sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa .```
```bwa index -a bwtsw gencode.v29.pc_transcripts.chr22.simplified.fa```

## Aligning the samples before shuffling to human chr 22 using bwa

```cd ~/workdir/Assignment/bwa_align```
for SAMPLE in 1 2 3 4 5;
    do
        R1=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_1.part_00${SAMPLE}.pe.trim.fastq
        R2=$HOME/workdir/Assignment/Mild_trimmed/SRR8797509_2.part_00${SAMPLE}.pe.trim.fastq
        /usr/bin/time -v bwa mem bwaIndex/gencode.v29.pc_transcripts.chr22.simplified.fa $R1 $R2 > SRR8797509_part_00${SAMPLE}.sam
done

## install the reqired tool to make alignment with hisat
##```conda install -c bioconda hisat2 ```

## Indexing
mkdir -p ~/workdir/Assignment/hisat_align/hisatIndex && cd ~/workdir/Assignment/hisat_align/hisatIndex
ln -s ~/workdir/sample_data/chr22_with_ERCC92.fa .

hisat2_extract_splice_sites.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > splicesites.tsv
hisat2_extract_exons.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92

## Aligning the samples after shuffling to human chr 22 using hisat 

```cd ~/workdir/Assignment/hisat_align```
for SAMPLE in 1 2 3 4 5;
    do
        R1=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_1.part_001.part_00${SAMPLE}.pe.trim.fastq
        R2=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_2.part_001.part_00${SAMPLE}.pe.trim.fastq
        hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuff_SRR8797509_part_00${SAMPLE}.sam
done


# Data Assembly 

## Assembly for the files resulted from BWA alignment on unshuffled data

```cd ~/workdir/Assignment/bwa_align```

## Prepare the SAM file for assembly
## install Samtools
## conda install samtools

## convert the SAM file into BAM file 
samtools view -bS SRR8797509_part_001.sam > SRR8797509_part_001.bam
samtools view -bS SRR8797509_part_002.sam > SRR8797509_part_002.bam
samtools view -bS SRR8797509_part_003.sam > SRR8797509_part_003.bam
samtools view -bS SRR8797509_part_004.sam > SRR8797509_part_004.bam
samtools view -bS SRR8797509_part_005.sam > SRR8797509_part_005.bam

## convert the BAM file to a sorted BAM file. 
samtools sort SRR8797509_part_001.bam -o SRR8797509_part_001.sorted.bam
samtools sort SRR8797509_part_002.bam -o SRR8797509_part_002.sorted.bam
samtools sort SRR8797509_part_003.bam -o SRR8797509_part_003.sorted.bam
samtools sort SRR8797509_part_004.bam -o SRR8797509_part_004.sorted.bam
samtools sort SRR8797509_part_005.bam -o SRR8797509_part_005.sorted.bam

## Export some useful statistics report for each sample indvidually
for f in 1 2 3 4 5;
    do
        samtools flagstat SRR8797509_part_00$f.sorted.bam > useful_stat_$f.txt;
done


## install required tool for assembly 
## install stringtie

## Assembly without known annotations
for SAMPLE in 1 2 3 4 5;
    do

        stringtie SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_free_${SAMPLE} -o ref_free_${SAMPLE}.gtf
done

# Assembly with known previous annotations
for SAMPLE in 1 2 3 4 5;
    do
        stringtie SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_sup_${SAMPLE} -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${SAMPLE}.gtf 
done

## Assembly for the files resulted from Hisat alignment on shuffled data

```cd ~/workdir/Assignment/hisat_align```

## convert the SAM file into BAM file 
samtools view -bS shuff_SRR8797509_part_001.sam > shuff_SRR8797509_part_001.bam
samtools view -bS shuff_SRR8797509_part_002.sam > shuff_SRR8797509_part_002.bam
samtools view -bS shuff_SRR8797509_part_003.sam > shuff_SRR8797509_part_003.bam
samtools view -bS shuff_SRR8797509_part_004.sam > shuff_SRR8797509_part_004.bam
samtools view -bS shuff_SRR8797509_part_005.sam > shuff_SRR8797509_part_005.bam

## convert the BAM file to a sorted BAM file. 
samtools sort shuff_SRR8797509_part_001.bam -o shuff_SRR8797509_part_001.sorted.bam
samtools sort shuff_SRR8797509_part_002.bam -o shuff_SRR8797509_part_002.sorted.bam
samtools sort shuff_SRR8797509_part_003.bam -o shuff_SRR8797509_part_003.sorted.bam
samtools sort shuff_SRR8797509_part_004.bam -o shuff_SRR8797509_part_004.sorted.bam
samtools sort shuff_SRR8797509_part_005.bam -o shuff_SRR8797509_part_005.sorted.bam

## Export some useful statistics report for each sample indvidually
for f in 1 2 3 4 5;
    do
        samtools flagstat shuff_SRR8797509_part_00$f.sorted.bam > shuff_useful_stat_$f.txt;
done


## install required tool for assembly 
## install stringtie

## Assembly without known annotations
for SAMPLE in 1 2 3 4 5;
    do

        stringtie shuff_SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_free_${SAMPLE} -o ref_free_${SAMPLE}.gtf
done

# Assembly with known previous annotations
for SAMPLE in 1 2 3 4 5;
    do
        stringtie shuff_SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_sup_${SAMPLE} -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${SAMPLE}.gtf 
done



# Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation.

##Create virtual evironment with conda

conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py


pypy3 -m pip install gffutils numpy tqdm 'intervaltree<3.0'


mkdir -p ~/workdir/Assignment/gtf-compare/gtfs && cd ~/workdir/Assignment/gtf-compare/gtfs
ln -s ~/workdir/Assignment/bwa_align/ref_free_*.gtf .
ln -s ~/workdir/Assignment/bwa_align/ref_sup_*.gtf .
mkdir -p ~/workdir/Assignment/gtf-compare/method_one && cd ~/workdir/Assignment/gtf-compare/method_one
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py

for f in 1 2 3 4 5;
        do
                pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/ref_free_*.gtf
                pypy3 stat.py
done
