# CNR Bioinformatics Workshop

[Back to curriculum](README.md)

# 3. Secondary Analysis
We are going to analyze aligned scRNAseq data in R, the most popular language for this type of work at the moment.

## 3.1 Install R libraries with Miniconda
First, let's install some packages with ```conda```. Fire up the cmd.
```
conda activate myEnv
conda install -c conda-forge r-base
conda install -c conda-forge r-seurat
conda install -c conda-forge r-ggplot2
conda install -c conda-forge r-dplyr
```

## 3.1 Import R libraries, set working directory and data
```
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

```


## 3.2 Quality control

```
df0[["percent.mt"]] <- PercentageFeatureSet(df0, features = intersect(mito$V1, rownames(df0)))
```
## 3.3 Principal component analysis

## 3.4 UMAP

## 3.5 Clustering analysis

## 3.6 Differential gene expression analysis


[<< Previous](2-alignment.md) [Next >>](4-visualization.md)  
