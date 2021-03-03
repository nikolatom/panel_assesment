if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("ggplot2", "RColorBrewer", "dplyr", "gridExtra", "maftools", "formattable"))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(maftools)
library(formattable)


############Plot overall coverage distribution#############
setwd("../results/bedtools/")
dir.create("../figures")
# get file names
print(files <- list.files(pattern="all.txt$"))

# load and parse the files
cov<-NULL
for (f in files){
  thisDF=read.table(f, sep="\t")[,2:5]
  colnames(thisDF)[c(1,4)]=c("coverage", "fraction")
  thisDF$fraction=1-cumsum(thisDF$fraction)
  cov =rbind(cov, cbind(thisDF, sample= gsub("^(SAMPLE[1-3]).*$", "\\1", f)))
}

# plot the data
maxCov=10000
p<-ggplot(data = subset(cov, coverage<maxCov), aes(x= coverage, y=fraction, colour=sample)) + ylim(0, 1)

plot1 = ggplot(data = subset(cov, coverage<maxCov), aes(x=coverage, y=fraction)) + 
  geom_line(aes(color = sample)) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold")) +
  ggtitle("Overall coverage depth and uniformity") +
  theme(plot.title = element_text(size = 20, face = "bold"))
#plot2 = ggplot(data = subset(cov, coverage>maxCov), aes(x=coverage, y=fraction)) + 
  #geom_line(aes(color = sample))
      
#grid.arrange(plot1, plot2, ncol=1)    

png(filename="../figures/coverage_fraction.png", width=1000, height=800)
print(plot1)
dev.off() 

############Print mean per target coverage#############

setwd("../picard_metrics")

list.files()

sample_S1 <- read.csv("AH_S1_L001.perTargetCov_amplicon.txt", header = TRUE, sep = "\t")
sample_S2 <- read.csv("CH_S2_L001.perTargetCov_amplicon.txt", header = TRUE, sep = "\t")
sample_S3 <- read.csv("CH_S2_L001.perTargetCov_amplicon_rmdup.txt", header = TRUE, sep = "\t")


sample_S1$target <- seq.int(nrow(sample_S1))

barplot(sample_S1$mean_coverage) 
cov_plot_S1<-ggplot(data=sample_S1, aes(y=mean_coverage, x=target)) +
  geom_bar(color = "violet", stat="identity")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold")) +
  ggtitle("Per target mean coverage AH_S1") +
  theme(plot.title = element_text(size = 20, face = "bold"))
png(filename="../figures/AH_S1_target_coverage.png", width=1000, height=800)
print(cov_plot_S1)
dev.off() 

sample_S2$target <- seq.int(nrow(sample_S2))
barplot(sample_S2$mean_coverage) 
cov_plot_S2<-ggplot(data=sample_S2, aes(y=mean_coverage, x=target)) +
  geom_bar(color = "green", stat="identity")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold")) +
  ggtitle("Per target mean coverage CH_S2") +
  theme(plot.title = element_text(size = 20, face = "bold"))
png(filename="../figures/CH_S2_target_coverage.png", width=1000, height=800)
print(cov_plot_S2)
dev.off() 

sample_S3$target <- seq.int(nrow(sample_S3))
barplot(sample_S3$mean_coverage) 
cov_plot_S3<-ggplot(data=sample_S3, aes(y=mean_coverage, x=target)) +
  geom_bar(color = "orange", stat="identity")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold")) +
  ggtitle("Per target mean coverage CH_S3 without duplicates") +
  theme(plot.title = element_text(size = 20, face = "bold"))
png(filename="../figures/CH_S3_target_coverage_rmdup.png", width=1000, height=800)
print(cov_plot_S3)
dev.off() 

png(filename="../figures/3panels_per_target.png", width=1000, height=800)
grid.arrange(cov_plot_S1, cov_plot_S2, cov_plot_S3, ncol=1)  
dev.off() 

mean(sample_S1$mean_coverage)
sd(sample_S1$mean_coverage)
min(sample_S1$mean_coverage)
max(sample_S1$mean_coverage)

mean(sample_S2$mean_coverage)
sd(sample_S2$mean_coverage)
min(sample_S2$mean_coverage)
max(sample_S2$mean_coverage)

mean(sample_S3$mean_coverage)
sd(sample_S3$mean_coverage)
min(sample_S3$mean_coverage)
max(sample_S3$mean_coverage)

################MAF TOOLS ANALYSYS#################

setwd("../varcall/maf")
#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = "allsamples.maf", clinicalData = laml.clin, removeDuplicatedVariants = FALSE, useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL, gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = "all", cnTable = NULL, isTCGA = FALSE, vc_nonSyn = c("3'Flank", "5'Flank","Silent","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation","Intron","3'UTR","5'UTR","Splice_Region"), verbose = TRUE)
vc_cols = brewer.pal(12,"Paired")

#Shows sample summry. ##NOT WORKING
df = getSampleSummary(laml)

##print table with maf summary

varTable<-formattable(df, 
                      align = c("l",rep("r", NCOL(df) - 1)),
                      list(`Tumor_Sample_Barcode` = formatter("span", style = ~ style(color = "grey", font.weight = "bold")), 
                           area(col = 2:10) ~ color_tile("#DeF7E9", "#71CA97")))

###save manually
#png(filename="../../figures/varTable.png")
print(varTable)
#dev.off() 


