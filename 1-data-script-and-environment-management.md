[Back to curriculum](README.md)
# CNR Bioinformatics Workshop

## 1. Package management
A few to many softwares need to be installed in order to accomplish any and every bioinformatics analysis task. If managed manually, it will quickly become messy and tedious, due to dependency and conflict issues between softwares. A software (aka library, aka package, aka module) management...software helps automating the process, and thus removing the headache.  
  
#### conda   
Conda is a light weight management software that is specifically suitable for biological studies.  
  
#### pip  
pip is another such software that may be used from time to time.  
   
### 1.1 Create a conda environment
Conda environments are independent "spaces" within a computer. Anything happens in one environment do not affect any other environment or outside of the conda "ecosystem", so it's a great way to keep the working environment tidy. Now, create your first conda environment with the following command:  
```
conda create -n myEnv
```
You can activate your environment with:  
```
conda activate myEnv
```
  
### 1.2 Install your first software in your new conda environment
You typically need to install various softwares within a conda environment to make it useful. Let's try with:
```
conda install -c conda-forge r-base
```

## 2. Softwares used to download data (on Mac or Linux)
#### wget  
[wget](https://www.gnu.org/software/wget/) is a software to download file (or files) from a link ([http](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol) or [ftp](https://en.wikipedia.org/wiki/File_Transfer_Protocol)).
For example, you can download Gencode mouse gene annotation v25 with the following code:
```
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
```
#### SRA Toolkit and fastq-dump/fasterq-dump  
**Sequence Read Archive (SRA)** is "the largest publicly available repository of high throughput sequencing data". **SRA Toolkit** are a suite of commandline softwares that alowing communication and interaction with the repository. **fastq-dump** is one of the commands within the toolkit used to download datasets from the database as fastq files.
**Fasterq-dump** is a stand-alone, modernized version of fastq-dump that allows multi-threading.

For the purpose of this workshop, we will use SRA Toolkit. Install it with conda:
```
conda install -c bioconda sra-tools
```

## 3. Code/Script management
1. Git
2. github

[<< Previous](0-the-basics.md)  [Next >>](2-alignment.md)  
