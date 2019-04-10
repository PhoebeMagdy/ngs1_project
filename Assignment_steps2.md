# Files Manipulation 

## 1. Download the files from SRA database
```wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra```

## 2. Download to the SRAtoolkit
```wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz```

## 3. unzip the SRAtoolkit
```tar -xzf sratoolkit.2.9.6-ubuntu64.tar.gz```

## 4. Add SRAtoolkit to the PATH:
```PATH=$PATH: /home/ngs-01/Downloads/sratoolkit.2.9.6-ubuntu64/bin```

## 5. Convert *.SRA file to FASTQ file using fastq-dump and downloading the first 5000000 reads only:
```fastq-dump - X 5000000 SRR8797509.sra```

## 6-counting number of sequences in the fastq file 
```grep '@' SRR8797509.fastq | wc -l```

## 7- separating each million in a separate file  
```cat SRR8797509.fastq | paste - - - - | awk 'NR>=1 && NR<=1000000' > sample1-1```

```cat SRR8797509.fastq | paste - - - - | awk 'NR>=1000001 && NR<=2000000' > sample2-1```

```cat SRR8797509.fastq | paste - - - - | awk 'NR>=2000001 && NR<=3000000' > sample3-1```

```cat SRR8797509.fastq | paste - - - - | awk 'NR>=3000001 && NR<=4000000' > sample4-1```

```cat SRR8797509.fastq | paste - - - - | awk 'NR>=4000001 && NR<=5000000' > sample5-1```

### or we can use the seqkit tool to split the main file into separate files
```seqkit split2  SRR8797509.fastq -s 1000000```

## 8- shuffleing the main fastq file containing the 5000000 reads
```seqkit shuffle SRR8797509.fastq > SRR8797509_shuffled.fastq```

## 9- resampling the shuffled file 
```cat SRR8797509_shuffled.fastq | paste - - - - | awk 'NR>=1 && NR<=1000000' > sample1-2```

```cat SRR8797509_shuffled.fastq | paste - - - - | awk 'NR>=1000001 && NR<=2000000' > sample2-2```

```cat SRR8797509_shuffled.fastq | paste - - - - | awk 'NR>=2000001 && NR<=3000000' > sample3-2```

```cat SRR8797509_shuffled.fastq | paste - - - - | awk 'NR>=3000001 && NR<=4000000' > sample4-2``` 

```cat SRR8797509_shuffled.fastq | paste - - - - | awk 'NR>=4000001 && NR<=5000000' > sample5-2 ```
## split the shuffled files 
```seqkit split2  SRR8797509_shuffled.fastq -s 1000000```

## checking the first sample before and after shuffling for quality control using FASTQC
```fastqc -t 1 -f fastq -noextract SSRR8797509.part_001.fastq SSRR8797509_shuffled.part_001.fastq```