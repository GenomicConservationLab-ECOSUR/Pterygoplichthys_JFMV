##### calculate SFS####
#####run SFS

library(vcfR)
library(dartR)
BiocManager::install("SNPRelate")
library(SNPRelate)

A1 <- read.vcfR("filtrovcf_maf0.5.recode.p.snps.vcf") #### for transcriptome data
gnlai <- vcfR2genlight (A1)
name <- read.table("nombres")
gnlai@pop <- as.factor(name$V2)
ploidy(gnlai) <- 2
gnlai

A1 <- read.vcfR("filtrovcf_maf0.5.recode.vcf") #### for de novo data
gnlai <- vcfR2genlight (A1)
name <- read.table("popmap")
gnlai@pop <- as.factor(name$V2)
ploidy(gnlai) <- 2
gnlai


gl.sfs(
  gnlai,
  minbinsize = 0,
  folded = TRUE,
  singlepop = FALSE,
  plot.out = TRUE,
  verbose = NULL
) -> sfs.pleco
