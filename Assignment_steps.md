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

## 8- shuffeling the main fastq file containing the 5000000 reads
```seqkit shuffle SRR8797509.fastq > SRR8797509_shuffeled.fastq```

## 9- resampling the shuffeled file 
```cat SRR8797509_shuffeled.fastq | paste - - - - | awk 'NR>=1 && NR<=1000000' > sample1-2```

```cat SRR8797509_shuffeled.fastq | paste - - - - | awk 'NR>=1000001 && NR<=2000000' > sample2-2```

```cat SRR8797509_shuffeled.fastq | paste - - - - | awk 'NR>=2000001 && NR<=3000000' > sample3-2```

```cat SRR8797509_shuffeled.fastq | paste - - - - | awk 'NR>=3000001 && NR<=4000000' > sample4-2``` 

```cat SRR8797509_shuffeled.fastq | paste - - - - | awk 'NR>=4000001 && NR<=5000000' > sample5-2 ```







