# Determining ancestry from SNP data 
# Author: Laura Budurlean

# Input: VCF files with SNV's called (using hg38) from your own samples
# Output: The final output will be a plot with your samples overlapped with super-populations from 1000 genomes data (https://www.internationalgenome.org/home) 

# You will need plink2, bcftools, vcftools

  # 1) bgzip and index sample VCF files then merge all into one vcf file.
srun bgzip file.vcf
srun bcftools index file.vcf.gz file2.vcf.gz
bcftools merge -m none -l vcfs.txt -Ov -o merged.vcf
  
  # 2) remove all duplicate mutation sites from merged.vcf
bcftools norm --rm-dup both -Oz -o merged_rm_dups.vcf.gz merged.vcf
  
  # 3) Convert merged_rm_dups.vcf.gz to plink format and also do some more VCF quality filtering. Output files with prefix "All" or another identifier of your choice. Depending on the plink version you are using if you receive a duplicate error  try increasing the default set-missing-var-id character limit from its default 23, using this --new-id-max-allele-len and set it to 100 or even higher until you dont get dups. This is a known issue in plink1.9 and plink2.0 solves it. Read more about it here: https://www.cog-genomics.org/plink/1.9/data#set_missing_var_ids
plink2 --vcf merged_rm_dups.vcf.gz --set-missing-var-ids @:#[b38]\$1\$2 --new-id-max-allele-len 250 --allow-extra-chr --chr 1-22 --double-id --max-alleles 2  --vcf-filter --vcf-min-gq 15 --make-bed --out All 

  # 4) Remove rare variants from "All" prefix files and create output files with prefix "All.geno" or another identifier of your choice. 
plink2 --bfile All --maf 0.1 --make-bed --out All.geno
# Note: Check if the All.geno.bim file has "chr" prefixes, if it does remove them: create a file with old names and new names and run plink2 with updated names list. You only need to do this if you see that you have "chr" prefixes in the All.geno.bim file. Otherwise skip this step and run plink normally without -update-name option.
#	awk '{print $2}' All.geno.bim > tmp1
#	sed 's/chr//g' tmp1 > tmp2
#	paste tmp1 tmp2 > updatenames.list
#	plink2 --bfile All.geno --update-name updatenames.list --snps-only --make-bed --out All.geno.snps

  # 5) Retain SNPs only, bgzip, and index. This file is what contains your samples variants that will be input into PCA later.
plink2 --bfile All.geno --snps-only --make-bed --recode vcf --out All.geno.snps
bgzip All.geno.snps.vcf
bcftools index All.geno.snps.vcf.gz  

  # 6) Create "site2.bed" file which is SNP positions present in your samples, that you will then use to filter out variants in the 1000 genomes samples. You need to make 2 files (1 for your samples and 1 for 1000 genomes) that only have the SNPs present in both files.
awk -F '\t' '{print $1"\t"$4-1"\t"$4}' All.geno.snps.bim > site2.bed

  # 7) Filter your sites against the 1000 genome sites. You will need phased shapeit2 files from 1000genomes "shapeit2_integrated..." (http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/)
# Set up an array as a bash script for every chromosome.
#SBATCH --array=1-22
array=(1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22)
number=${array[$SLURM_ARRAY_TASK_ID-1]}

bcftools filter -T site2.bed -Oz -o chr${number}.1kg.sites.vcf.gz ALL.chr${number}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 

  # 8) Combine all the chromosome files using vcftools, bgzip and index.
vcf-concat chr*.1kg.sites.vcf.gz > merged.1kg.sites.vcf
bgzip merged.1kg.sites.vcf
bcftools index merged.1kg.sites.vcf.gz

  # 9) Intersect the All.geno.snps.vcf.gz and merged.1kg.sites.vcf.gz. Set the # of threads as needed
bcftools isec --threads ## -p isec -n =2 All.geno.snps.vcf.gz merged.1kg.sites.vcf.gz
bgzip --threads 64 0000.vcf
bgzip --threads 64 0001.vcf 
bcftools index 0000.vcf.gz
bcftools index 0001.vcf.gz
bcftools merge -m none 0000.vcf.gz 0001.vcf.gz -Oz -o merged_ALL_1kg_samples.vcf.gz

  # 10) Set missing genotypes to 0|0
bcftools +setGT merged_ALL_1kg_samples.vcf.gz -- -t . -n 0p\n > merged_ALL_1kg_samples_setGT.vcf.gz
plink2 --vcf merged_ALL_1kg_samples_setGT.vcf.gz --set-missing-var-ids @:#[b38]\$r\$a --make-bed --out unpruned1kg_setGT
plink2 --bfile unpruned1kg_setGT --indep-pairwise 100kb 1 .1
plink2 --bfile unpruned1kg_setGT --extract plink2.prune.in --maf 0.1 --pca biallelic-var-wts --make-bed --out All.geno.prune.setGT
 
 # 11) PCA scoring. You may score as many principal components as needed. I recommend at least the first three for samples that map ambiguously when graphing PC1 and PC2. You may wish to graph PC1 and PC3 as well. Keeping in mind that the first two principal components will explain the most amount of variation in the data.
plink2 --bfile ${filedir}/unpruned1kg_setGT --score ${filedir}/All.geno.prune.setGT.eigenvec.var 2 3 5 --out 1kg-PC1-score
plink2 --bfile ${filedir}/unpruned1kg_setGT --score ${filedir}/All.geno.prune.setGT.eigenvec.var 2 3 6 --out 1kg-PC2-score
plink2 --bfile ${filedir}/unpruned1kg_setGT --score ${filedir}/All.geno.prune.setGT.eigenvec.var 2 3 7 --out 1kg-PC3-score

# We continue in R to plot our samples with 1000 genomes super-population data. See "plot.R" script
