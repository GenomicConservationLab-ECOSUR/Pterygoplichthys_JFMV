               ##### BASIC STATISTICS#####
# We load and read the vcf file and convert it to different formats
# activate the following libraries
library(adegenet)
library(vcfR)
library(hierfstat)
library(ape)
library(corrplot)

A1 <- read.vcfR("refmap.recode.p.snps.vcf")
genynd <- vcfR2genind(A1)# we convert to genind object normal analysis
pgn1 <- vcfR2genind(A1) ## to run exploratory analysis with the specimens (1,2 and 3)
gnlai <- vcfR2genlight (A1)
name <- read.table("popmap")
genynd@pop <- as.factor(name$V2) #we establish the name of the town
genynd
gnlai@pop <- as.factor(name$V2)
ploidy(gnlai) <- 2
gnlai
# Analisis basicos de diversidad genetica #
###########################################

divers <- basic.stats(genynd) # Global analysis, all populations
divers$overall # Average of all data. All diversity measures

# To obtain the averages per population
colnames(divers$Ho) # To see the name of populations
levels(pgn@pop) 

nombresPobs <- colnames(divers$Ho) # Make an object with the name of the populations
HetObs <-  colMeans(divers$Ho)   # Average observed heterozygosity of the population 
Richness <- allelic.richness(Genind1) # Estimate allelic richness (microsatellites), rarefied allelic counts, by locus and population

# Make a bar chart with the observed heterozygosity values.
newnomb <- c("Catazajá", "Grijalva", "Lacantún", "Tres Brazos", "Tenosique")
ggplot(HetObs, xlab = "Poblaciones", ylab="Heterosigocidad observada (HO)", col="gray", names.arg = newnomb)
abline (h=0.2357) # If we want to add a horizontal line that shows the average value, taking the mean of the divers$overall table

##########################
# We can do the same to obtain the average values of the other diversity measures.
##########################

# dev.new()

par(mfrow=c(1,3)) 
myCol <- c("blue","yellow","green","grey","red")
barplot(colMeans(divers$Hs, na.rm = T), col = myCol , names.arg = newnomb, las=2, axis.lty = 6, axisnames = T, ylab = "H.S.", xlab = "Populations")
barplot(colMeans(divers$Ho, na.rm = T), col = myCol, names.arg = newnomb, las=2, axis.lty = 6, axisnames = T, ylab = "H.O.",  xlab = "Populations")
barplot(colSums(Richness$Ar),  las=2, col = myCol, ylab = "N. Alleles")

###################################
# Genetic structure analysis #
###################################

# We can estimate distances to heatmap genetic distance between individuals.
x.dist <- dist(pgn1)
heatmap(as.matrix(x.dist))

# Fst paired between populations
a.fstat <- genind2hierfstat(b)
Fstat.metrics <- pairwise.neifst(a.fstat)
par(mfrow=c(1,1)) 
corrplot(as.matrix(Fstat.metrics), is.corr = F, type = "lower", diag=F) 

# Estimate Nei's genetic distance between populations
pgp <- genind2genpop(pgn1) # We first need to transform our data to genpop format.
y.dist<- dist.genpop(pgp, method = 1)
heatmap(as.matrix(y.dist), symm = T)
min(y.dist, na.rm = T)
max(y.dist, na.rm = T)

# Make a UPGMA of how similar the populations are to each other.
c <- hclust(y.dist, method = "average")
par(mfrow=c(1,2))
plot(c, main= "UPGMA dendrogram", xlab = "Pterygoplichthys spp")
corrplot(as.matrix(Fstat.metrics), is.corr = F, type = "lower", diag=F, method = "square", col.lim = c(0,1))
# Make a DAPC with colors to identify genetic groups. This can be done with the populations already defined
par(mfrow=c(1,1))
dpca1 <- dapc() #I retained 2 discriminants
dpca3 <- dapc(a)
dpca4 <- dapc(a)

scatter.dapc(dpca1, cex=4, scree.da =FALSE)
compoplot.dapc(dpca1)

# calculate frequency spectrum

myfreg <- glMean(gl1)
hist(myfreg, proba=TRUE)

myfreg2 <- c(myfreg, 1-myfreg)
hist(myfreg2, proba=TRUE)

# Now we convert the object (A1) to genind and genlight format.

Genind1 <- vcfR2genind(A1)
gl1 <- vcfR2genlight(A1)

class(gl1)
str(gl1)
# explorar la base
gl1@ind.names
popmap <- read.table("popmap")
genind1@pop <- as.factor(popmap$V2) #we assign populations to the data
genlai@pop <- as.factor(popmap$V2)
gl1@pop
Genind1@pop

genpop1 <- genind2genpop(genind1) # we transform to genpop format
pgnpop <- genind2genpop(pgn)
class(Genind1_pop)
Genind1_pop
# cargamos librerías
library(adegenet)
library(pegas)
library(hierfstat)
########################################################################
glPlot(gl12, posi="topleft")
gl12$ind.names
gl12$loc.names[1:10]
glPlot(gl1, posi="bottomleft")
pca1 <- glPca(denovo_gl, nf=4)
pca1 
plot(pca1$scores[,1], pca1$scores[,2], #plot individuals
     cex=2, pch=20, col=denovo_gl$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 3095 SNPs)") 
legend("topleft", 
       legend=unique(denovo_gl$pop), 
       pch=20)
### testing basic statistics
summary(Genind1.pop)
div <- summary(genpop1)
names(div)
### with hierfstats
pgn.hfstat <- genind2hierfstat(pruenagn)
basicstat <- basic.stats(pgn.hfstat, diploid = TRUE, digits = 2) 
names(basicstat)
boot.ppfis(pgn.hfstat)
x <- indpca(pgn.hfstat) 
plot(x, cex = 0.7)
library(pegas)
hw.test(pgn, B = 1000)
## see explained variance
var_frag <- pca1$eig/sum(pca1$eig)
signif(sum(var_frag[1:2])*100,2)
#dapc
pnw.dapc <- dapc(gl12)
# plot result
scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)
compoplot(pnw.dapc, posi = 'topright')
loadingplot(pnw.dapc$loadings)
loadingplot(pca1)
#mejoramos el complot
library(tidyr)
#change name of population
gl12 <- gl1
gl12@pop
gl12@pop <- as.factor(name.loc$V2) #assign names
Genind1@pop <- as.factor(name.loc$V2)
# we convert file .genind to .genepop
Genind1.pop <- genind2genpop(Genind1)
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl12)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
# volvemos a graficar
p <- ggplot(dapc.results, aes(x=indNames, y=value, fill=pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(aesthetics = "colour") 
p <- p + facet_grid(~pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

### view statistics by population
library(mmod)
name.loc <- read.table("popmap.txt")
Genind12 <- Genind1
strata(Genind12) <- name.loc
Genind12
setPop(Genind12) <- ~V2 # Use @pop to analyze by population
diff_pop <- diff_stats(Genind12)
diff_pop
## we plot results
library(ggplot2)
library(reshape2)
per.locus <- melt(diff_pop$per.locus, varnames = c("Locus", "Statistic"))
stats     <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
glob      <- data.frame(Statistic = stats, value = diff_pop$global)
head(per.locus)
head(glob)
# graficamos
ggplot(per.locus, aes(x = Statistic, y = value)) +
  geom_boxplot() +
  geom_point() +
  geom_point(size = rel(3), color = "red", data = glob) +
  ggtitle("Estimates of population differentiation")
per.locus

### explore results
basicstat$overall
basicstat$n.ind.samp
basicstat$Fis
basicstat$perloc
## arbol de distancia
x.dist <- dist(x)
x.dist <- poppr::bitwise.dist(x)
tree <- aboot(gl12, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

# coloring the tips of the tree
cols <- brewer.pal(n = nPop(gl12), name = "Paired")
plot.phylo(tree, cex = 0.8, font = 1, adj = 0, tip.color =  cols[pop(gl12)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 1, xpd = TRUE)
#legend(35,10,c("Cat","Grij","Lac"...),cols, border = FALSE, bty = "n")
legend('topright', legend = c("Cat","Grij","Lac", "TB", "Teno"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#########################################################################################

myFreq <- glMean(gl12)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=2)


tre <- nj(dist(as.matrix(gl12)))
tre
plot(tre, typ="fan", cex=0.7)
title("NJ tree of the US influenza data")
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of the US influenza data")

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

#### we use the bartlet test for global HW

names(divers.end)
bartlett.test(list(div$Hexp, div$Hobs))