[Back to curriculum](README.md)  

## 1. Cellranger
When dealing with single cell RNA- or ATAC-seq generated with the 10X Genomics platform, Cellranger is recommend because of its simplicity. The Cellranger program only runs on Linux system. For Mac and PC users, we will utilize the [cloud analysis](https://www.10xgenomics.com/products/cloud-analysis) platform.

## 2. STAR aligner  
One of the first decisions one has to make when it comes to alignment (aka "mapping") is the choice of aligner. There are more aligners than anyone can care about out there but I personally use STAR when dealing with RNAseq data (sometimes even single cell RNAseq data). For any alignment task, regardless of the aligner, we need <b> a genome reference </b> and <b> an annotation file </b> of the corresponding species that we want to align to. You can download these from various sources (Ensembl, UCSC, Genconde etc.).  As an example, here is some code I used to download Human genome reference release 38 (GRCh38.p13) and Gencode gene annotation, which can be find on [this page](https://www.gencodegenes.org/human/)
```
wget -o /path/to/genome_references/GRCh38 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget -o /path/to/genome_references/GRCh38 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
```
Once downloaded, uncompress them with ```gunzip``` and they are ready to be used.
```
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v38.annotation.gtf.gz
```

Next, STAR requires a pre-built index to run, which is built using the ```genomeGenerate``` run mode:
```
STAR --genomeDir /path/to/genome_references/GRCh38 \
--runMode genomeGenerate \
--genomeFastaFiles /path/to/genome_references/GRCh38/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /path/to/genome_references/GRCh38/gencode.v38.annotation.gtf \
--runThreadN 10
```
 
The out put of the above program is a whole bunch of indexes that we don't need to worry about, but are required in the next step.

Finally, we can run the alignment step with:
```
for file in *.fastq.gz; do

prefix=${file%%.fastq.gz}

STAR --genomeDir /path/to/genome_references/GRCh38 \
     --runThreadN 10 \
     --genomeLoad LoadAndKeep \
     --readFilesIn $file \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 30000000000 \
     --outTmpDir ~/${prefix}_tmp \
     --outFileNamePrefix /path/to/out_put_dir

done

STAR --genomeDir /path/to/genome_references/GRCh38 --genomeLoad Remove 
```
## 3. (Optional) Quality control with fastqc
Download [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) and complete installation following instructions.

## 4. 10X Genomics cloud analysis
Register an account at: https://cloud.10xgenomics.com/signin
