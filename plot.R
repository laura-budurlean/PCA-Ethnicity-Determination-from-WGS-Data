# Load libraries
library(ggplot2)
library(ggrepel)

# Reading in the sample info file. This file contains your own sample information and information on 1000 genomes individuals. Your samples should be added after all the 1000 genomes samples in order to show up as the top layer when plotting. The file should be CSV formatted as so: 
# Sample,Family ID,Population,Population Description,Gender,Relationship,Unexpected Parent/Child ,Non Paternity,Siblings,Grandparents,Avuncular,Half Siblings,Unknown Second Order,Third Order,Other Comments
# HG00096,HG00096,GBR,British in England and Scotland,male,,,,,,,,,,
# your_sample_ID,your_sample_ID,your_cohort_name,,,,,,,,,,,,

kg<-read.csv("C:/.../.../20130606_sample_info.csv",header=T,stringsAsFactors=F)

# Creating another column at the end of the sample_info.csv file called "SuperPopulation" which is going to have 3 letter abbreviations of super-populations from the sub-populations. These super=pops are: EAS, EUR, AFR, AMR, SAS 
kg$SuperPopulation[kg$Population=="CHB"|kg$Population=="JPT"|kg$Population=="CHS"|kg$Population=="CDX"|kg$Population=="KHV"]<-"EAS"
kg$SuperPopulation[kg$Population=="CEU"|kg$Population=="TSI"|kg$Population=="FIN"|kg$Population=="GBR"|kg$Population=="IBS"]<-"EUR"
kg$SuperPopulation[kg$Population=="YRI"|kg$Population=="ASW"|kg$Population=="ACB"|kg$Population=="LWK"|kg$Population=="GWD"|kg$Population=="MSL"|kg$Population=="ESN"]<-"AFR"
kg$SuperPopulation[kg$Population=="MXL"|kg$Population=="PUR"|kg$Population=="CLM"|kg$Population=="PEL"]<-"AMR"
kg$SuperPopulation[kg$Population=="GIH"|kg$Population=="PJL"|kg$Population=="BEB"|kg$Population=="STU"|kg$Population=="ITU"]<-"SAS"
kg$SuperPopulation[kg$Population=="your_cohort_name"]<-"your_cohort_name"

# Reading in all the PCA score profile files for 1000 genomes data. Check headers and adjust headers= option as needed.
pc1<-read.table("1kg-PC1-score.sscore",header=T,stringsAsFactors=F)
pc2<-read.table("1kg-PC2-score.sscore",header=T,stringsAsFactors=F)
pc3<-read.table("1kg-PC3-score.sscore",header=T,stringsAsFactors=F)

# Create a data frame with the required information for plotting
kg2<-data.frame(FID=pc1$FID,IID=pc1$IID,PC1=pc1$SCORE1_AVG,PC2=pc2$SCORE1_AVG,PC3=pc3$SCORE1_AVG,SuperPopulation=kg$SuperPopulation)

# Plot your samples with 1000 genomes data. This is an example plot with labels. You may need to adjust these options as needed until you obtain your desired plot formats. 
ggplot(kg2,aes(x=PC1,y=PC2,color=SuperPopulation,shape=SuperPopulation, label=IID)) + geom_point(size = 3, shape = 16) + geom_text_repel(data = subset(kg2, SuperPopulation == "cohort"), size= 3, colour="black", max.overlaps = 100) + scale_color_manual(values=c("purple","orange","black","green","blue","red"))

