## Overview
This is a curriculum of a practical introduction workshop to bioinformatics analysis.  
This workshop is an very general and very basic introduction to tools and practices of bioinformatics analyses, focusing on RNAseq and single cell RNAseq. It is designed for an audience with limited to no bioinformatics or coding experience and hopefully will help them overcome the innate fear of the command line. Since all the contents as of now are derived exclusively from my own personal experience, I cannot claim them to be the most common, most efficient or best practices, but they got the jobs done for me when I, a neuroscientist with no formal training, needed. Once again, the goal isn't to provide a rigorous guide of standard practices but rather helping those who are in need to kick the can down the road.  

## Prerequisite
Miniconda: https://docs.conda.io/en/latest/miniconda.html

#### Miniconda installation  
Refer to the [installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) page for detail.
   * On Windows:  
   Download installer from [here](https://docs.conda.io/en/latest/miniconda.html)
   * On Mac or Linux:  
   Check [here](https://docs.conda.io/en/latest/miniconda.html)   
   Download installer with wget and install with bash:    

Mac
   ```
   wget Miniconda3-latest-MacOSX-x86_64.sh
   bash Miniconda3-latest-MacOSX-x86_64.sh
   ```
   Linux
   ```
   wget Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ``` 
   Follow prompts and finish the installation

## Topics
[0. The basics](0-the-basics.md)  
1. The Mighty CMD: command line interface  
2. Version control  

[1. Data, script and environment management](1-data-script-and-environment-management.md)
1. wget
2. fastq-dump
3. conda
4. Git/github

[2. Alignment and read quantification](2-alignment.md)
1. cellranger
2. The STAR of the show: Alignment with STAR/featureCounts

[3. Secondary analysis](3-secondary-analysis.md)
1. Getting started with R
2. Load libraries and datasets
3. Dimension reduction
4. Clustering
5. Differential gene expression

[4. Visualization](4-visualization.md)
1. ggplot2
2. Cerebro

## Not covered but useful
- Access high performance computer with PuTTy (on windows) or ssh (on Mac and Linux);
- Upload and download files with Filezilla;
- Upload and download files with scp;
- Parallelization: run a program with more than one processor.

