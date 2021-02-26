# PCA-Ethnicity-Determination-from-WGS-Data
A pipeline utilizing TCGA data and WGS data from your own samples to determine or validate ethnicity of an individual.

The goal of this pipeline is to determine ethnicity of an individual using sequencing data (SNPs) starting with variant called files (VCF) from those individuals. The cohort data is then combined/overlayed with TCGA data and PCA analysis is performed. PCA scores are then plotted along with TCGA data to provide a visual representation of where each individual falls on the overall PCA plot of ethnicities.

Some requirements for this pipeline.
- filtered VCF files 
- bcftools
- plink
- vcftools
- R and R-studio (for plotting)
- if processing many samples: high performance computing cluster
