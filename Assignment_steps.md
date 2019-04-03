# Files Manipulation 

## 1. Download the files from SRA database
```wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra```

## 2. Download to the SRAtoolkit
```wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz```

## 3. unzip the SRAtoolkit
```tar -xf sratoolkit.2.9.6-ubuntu64.tar.gz```

## 4. Add SRAtoolkit to the PATH:
```PATH=$PATH: /home/ngs-01/Downloads/sratoolkit.2.9.6-ubuntu64/bin```

## 5. Convert *.SRA file to FASTQ file using fastq-dump:
```fastq-dump SRR8797509.sra```



