install.packages("devtools") # First we install
library("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")

library("pcadapt")
library("qvalue")
library("OutFLANK")
library("ggplot2")

####### Part 1: identifying atypical snp with pcadap #### 
library(LEA)
library(lfmm)
library(qvalue)

path_to_file <- "D:/PLECO/Transcriptoma/Seleccion/LFMM/data.lfmm" #### we load the base in vcf or lfmm format
filename <- read.pcadapt(path_to_file, type = "lfmm")

x <- pcadapt(input = filename, K = 2) ###we make a pc
plot(x, option = "screeplot")  ## we plot it and see where it explains the largest variation

### we color the pc by location
plot(x,option="scores",pop=meta$Pop)

### now we will see which of the SNPs are driving the overall variation.
plot(x,option="manhattan")

#### see SNPs associated with a single component
x_cw <- pcadapt(filename,K=2,method ="componentwise")
plot(x_cw,option="manhattan",K=2)
snp_unpc <- get.pc(x_cw, outliers)
head(snp_unpc)

#### now to see statistical outliers we adjust for p
qval <- qvalue(x$pvalues)$qvalues
outliers <- which(qval<0.05)
length(outliers)

#### in the adjustment we will use the Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.05
outliers_bonferroni <- which(padj < alpha)
length(outliers_bonferroni)

#####obtain the list of atypical SNPs
snp_pc <- get.pc(x, outliers_bonferroni)
head(snp_pc)

meta <-read.table("meta.csv",h=T,sep = ",")

################now we will see graphically if there are atypical SNPs associated with salinity and RQI.
plot(x$scores[,1]~meta$Sal,pch=19,col="gray") # sal pc1
plot(x$scores[,2]~meta$Sal,pch=19,col="gray") # sal pc2
plot(x$scores[,1]~meta$RQI,pch=19,col="gray") # RQI pc1
plot(x$scores[,2]~meta$RQI,pch=19,col="gray") # RQI pc2

################################################################################
################################################################################
########## We now use Fst to identify outlayers#########################
################################################################################
library(OutFLANK)
library(vcfR)
base<-read.vcfR("D:/PLECO/Transcriptoma/Seleccion/LFMM/data.vcf")
gen <- extract.gt(base)
dim(gen)
head(gen[,1:10])
g<-gen
#####we convert to OutFlank format
g[gen %in% c("0/0")] <- 0
g[gen  %in% c("0/1")] <- 1
g[gen %in% c("1/1")] <- 2
g[is.na(g)] <- 9
tg <- t(g)
dim(tg)
levels(factor(as.numeric(unlist(tg))))
meta<-read.table("meta1.csv",h=T,sep = ",")
subpops <- c("GRIJ","TEN","LAC","CAT","TB")
subgen <- tg[meta$Pop%in%subpops,] #subset method 1

submeta <- subset(meta,Pop%in%subpops) #subset method 

identical(rownames(subgen),as.character(submeta$ID))

fst <- MakeDiploidFSTMat(subgen,locusNames=1:ncol(subgen),popNames=submeta$ID)



head(fst)
hist(fst$FST,breaks=50)
summary(fst$FST)
############ahora vemos valores atipicos estadisticos
OF <- OutFLANK(fst,LeftTrimFraction=0.01,RightTrimFraction=0.01,
               Hmin=0.05,NumberOfSamples=2,qthreshold=0.01)
OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.05,titletext=NULL)

P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.05,Hmin=0.1)
outliers <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers)
######we plot the atypical snp with a Manhatan
plot(P1$LocusName,P1$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1$LocusName[outliers],P1$FST[outliers],col="magenta")


################ Part 2: we run with LFMM ###########

library(LEA)

vcf2lfmm("D:/PLECO/Transcriptoma/Seleccion/LFMM/data.vcf", "data_lfmm") # Transform database to lfmm format
proyecto <- lfmm("data.lfmm", "amb.env", K=2, repetitions=10, 
                 all=TRUE, missing.data=TRUE,iterations=50000, burnin=5000, random.init=TRUE, 
                 project="new") # Correrlo desde R consola
#The project is saved into :
  # data_amb.lfmmProject 

#To load the project, use:
  project = load.lfmmProject("data_amb.lfmmProject")

#To remove the project, use:
  remove.lfmmProject("data_amb.lfmmProject")

#[1] "*************************************"
#[1] "* lfmm K = 2  repetition 1  all     *"
#[1] "*************************************"
#Summary of the options:
  
# -n (number of individuals)      103
# -L (number of loci)             3059
# -K (number of latent factors)   2
# -o (output file)                data_amb.lfmm/K2/run1/data_r1
# -i (number of iterations)       50000
# -b (burnin)                     5000
# -s (seed random init)           2129367179
# -p (number of processes (CPU))  1
# -x (genotype file)              data.lfmm
# -v (variable file)              amb.env
# -D (number of covariables)      3
# -d (the dth covariable)         1
# -a (all variable at the same time)
# -m (missing data)                 

  p1 <- lfmm.pvalues(proyecto, K = 2, d=2)
  p2 <- lfmm.pvalues(proyecto, K = 2, d=3)
###we calculate candidate spn as a function of FDR
  
  library(qvalue)
  
  alpha <- 0.05
  
  qval1 <- qvalue(p1$pvalues)$qvalues
  candidates1 <- which (qval1 < alpha)
  length (candidates1)
  
  qval2 <- qvalue(p2$pvalues)$qvalues
  candidates2 <- which (qval2 < alpha)
  length (candidates2)
  ###we calculate with bonferroni correction
  padj1 <- p.adjust(p1$pvalues,method="bonferroni")
  alpha <- 0.05
  outliers_sal <- which(padj1 < alpha)
  length(outliers_sal)
  
  padj2 <- p.adjust(p2$pvalues,method="bonferroni")
  alpha <- 0.05
  outliers_rqi <- which(padj2 < alpha)
  length(outliers_rqi)
  
    # Make Manhattan plots for each variable highlighting candidates.

  plot (-log10(p1$pvalues), type="p", col="grey", cex=.4, pch=19, xlab="SNPs", ylab="-log10(p-value)", main="Sal")
  points(candidates1, -log10(p1$pvalues[candidates1]), col = "red")
  plot (-log10(p2$pvalues), col="grey", cex=.4, pch=19, xlab="SNPs", ylab="-log10(p-value)", main="RQI")
  points(candidates2, -log10(p2$pvalues[candidates2]), col = "red")
  
    # Identify the name of these loci
  d <- read.table ("posicionesNuevas.txt", header=FALSE)
  
  Posi1 <- d[candidates1, ]
  Posi2 <- d[candidates2, ]

  
  