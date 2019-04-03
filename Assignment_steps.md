# Files Manipulation 

## Download the files from SRA database
```wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra```

## Download to the SRAtoolkit
```wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz```

## unzip the SRAtoolkit
```tar -xf sratoolkit.2.9.6-ubuntu64.tar.gz```

## Add SRAtoolkit to the PATH:
```PATH=$PATH: /home/ngs-01/Downloads/sratoolkit.2.9.6-ubuntu64/bin```

## Convert *.SRA file to FASTQ file using fastq-dump:
```fastq-dump SRR8797509.sra```



