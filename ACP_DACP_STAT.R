#### we load the following libraries
library(adegenet)
library(vcfR)
library(hierfstat)
library(ape)
library(corrplot)
library(ggplot2)
library(RColorBrewer)


a1 <- read.vcfR("filtrovcf_maf0.5.recode.p.snps.vcf")
denovo_gl <- vcfR2genlight(a1)
ploidy(denovo_gl) <- 2
denovo_gnd <- vcfR2genind(a1) #convert to genind
name.loc <- read.table("nombres")
denovo_gl@pop <- as.factor(name.loc$V2)
denovo_gnd@pop <- as.factor(name.loc$V2)
denovo.pca <- glPca(c, nf = 3)
par(mfrow = c(1,1))
barplot(100*denovo.pca$eig/sum(denovo.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

denovo.pca.scores <- as.data.frame(denovo.pca$scores)
denovo.pca.scores$pop <- pop(c)

transc.pca <- glPca(c, nf = 3)
barplot(100*transc.pca$eig/sum(transc.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

transc.pca.scores <- as.data.frame(transc.pca$scores)
transc.pca.scores$pop <- pop(c)

library(ggplot2)

set.seed(9)
par(mfrow = c(1,2))
p <- ggplot(denovo.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = c("#0000FF", "#458B00", "#030303", "#FFD700","#8B0000")) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p

set.seed(9)
p1 <- ggplot(transc.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p1 <- p1 + geom_point(size=2)
p1 <- p1 + stat_ellipse(level = 0.95, size = 1)
p1 <- p1 + scale_color_manual(values = c("#0000FF", "#458B00", "#030303", "#FFD700","#8B0000")) 
p1 <- p1 + geom_hline(yintercept = 0) 
p1 <- p1 + geom_vline(xintercept = 0) 
p1 <- p1 + theme_bw()

p1

cols <- c("#0000FF", "#458B00", "#030303", "#FFD700","#8B0000")
denovo.dapc <- dapc(c, n.pca = 3, n.da = 2)
scatter(denovo.dapc, col = cols , cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)

transc.dapc <- dapc(cgl, n.pca = 3, n.da = 2)
scatter(transc.dapc, col = cols , cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)


compoplot(denovo.dapc,col = cols, posi = 'top')
compoplot(transc.dapc,col = cols, posi = 'top')

library(tidyr)
#denovo
dapc_denovo.results <- as.data.frame(denovo.dapc$posterior)
dapc_denovo.results$pop <- pop(c)
dapc_denovo.results$indNames <- rownames(dapc_denovo.results)
dapc_denovo.results <- pivot_longer(dapc_denovo.results, -c(pop, indNames))
head(dapc_denovo.results, n= 6)
colnames(dapc_denovo.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
pdn <- ggplot(dapc_denovo.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
pdn <- pdn + geom_bar(stat='identity') 
pdn <- pdn + scale_fill_manual(values = c("#0000FF", "#458B00", "#030303", "#FFD700","#8B0000")) 
pdn <- pdn + facet_grid(~Original_Pop, scales = "free")
pdn <- pdn + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
pdn

#transcriptoma
dapc_transc.results <- as.data.frame(transc.dapc$posterior)
dapc_transc.results$pop <- pop(cgl)
dapc_transc.results$indNames <- rownames(dapc_transc.results)
dapc_transc.results <- pivot_longer(dapc_transc.results, -c(pop, indNames))
head(dapc_transc.results, n= 6)
colnames(dapc_transc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

ptp <- ggplot(dapc_transc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
ptp <- ptp + geom_bar(stat='identity') 
ptp <- ptp + scale_fill_manual(values = c("#0000FF", "#458B00", "#030303", "#FFD700","#8B0000")) 
ptp <- ptp + facet_grid(~Original_Pop, scales = "free")
ptp <- ptp + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
ptp

####calculate fst
#convert to vcfr to genind

denovo_genind <- vcfR2genind(a)
transc_genind <- vcfR2genind(a1)

###we recalculate statistics to run a comparison analysis

denovo_gnd.hfstat <- genind2hierfstat(denovo_gnd)

basicstat <- basic.stats(denovo_gnd, diploid = TRUE, digits = 2) 
names(basicstat)
boot.ppfis(Genind1.hfstat)