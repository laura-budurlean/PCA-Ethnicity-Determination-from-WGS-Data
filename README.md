# PCA-Ethnicity-Determination-from-WGS-Data

A pipeline utilizing TCGA data and WGS data from your own samples to determine or validate ethnicity of an individual.

The goal of this pipeline is to determine ancestry of an individual using sequencing data (SNPs) starting with hg38 variant called files (VCF) from those individuals. The cohort data is then combined/overlayed with 1000 genomes data and PCA analysis is performed. PCA scores are then plotted along with 1000 genomes data to provide a visual representation of where each individual falls on the overall PCA plot of ancestry. 

Some requirements for this pipeline.
- filtered VCF files from your own samples (hg38)
- bcftools
- plink2
- vcftools
- R and R-studio (for plotting)
- 1000 genomes data files http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/
- if processing many samples: high performance computing cluster

The output of this ancestry calling pipeline will give you a plot with 1000 genomes super populations and your own samples overlayed on top of the super population they most closely resemble based on the SNV data.

![example_PCA_for_github](https://github.com/laura-budurlean/PCA-Ethnicity-Determination-from-WGS-Data/assets/30268603/f57e7e1f-81b9-4f2e-b6e2-46e4e79a6f82)
