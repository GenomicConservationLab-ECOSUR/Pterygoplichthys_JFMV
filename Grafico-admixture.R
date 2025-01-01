####we plot values of k
K1 <-read.table(file = "Pleco_Reducido_Kerror_trn.txt", header = TRUE)
head(Q2)
plot(K1$X0.58149, data$valor, type= "b", col= "darkorchid3")

## grafico k
k10 <- K1[2,]
k10
K <- K1[-2,]
K
K <- rbind(K, k10)

plot(K$valor, xlab ="K value", ylab = "CV error")

library(ggplot2)
ggplot(a$V4, aes(x=K, y=valor))+
  geom_line(colour= "black") +
  geom_point(size=2, shape= 21, fill="white", colour= "red") +
  theme_minimal()

plot(K, pch= 19, xlab= "K value", ylab= "CV error")

### plot data we call from folder###########
Q2 <-read.table(file = "nm.plink.2.Q", header = FALSE)
Q3 <-read.table(file = "nm.plink.3.Q", header = FALSE)
### graficar datos
colnames(Q2) <- c("Cluster 1", "Cluster 2")
Q2
colnames(Q3) <- c("Cluster 1", "Cluster 2", "Cluster 3")
DatPob <- read.table("popmap.txt", sep = "\t")
DatPob <- read.table("popmap1.txt", sep = "\t")
DatPob <- read.table("popmap1.txt", header = FALSE)
DatPob <- read.table(file= "popmap1.txt", header = FALSE)
DatPob <- read.table(file = "popmap1.txt", header = FALSE)
####graficamos valores de k
DatPob <-read.table(file = "popmap1", header = FALSE)
DatPob
colnames(DatPob) <- c("Ind", "Pop")
colnames(DatPob)
Q2 <- data.frame(DatPob$Ind, DatPob$Pop, Q2)
head (Q2)
colnames(Q2) <- c("Ind", "Pop", "Cluster 1", "Cluster 2")
head(Q2)
Q3 <- data.frame(DatPob$Ind, DatPob$Pop, Q3)
head (Q3)
colnames(Q3) <- c("Ind", "Pop", "Cluster 1", "Cluster 2", "Cluster 3")
head(Q3)
mdat3 <- melt (Q3, id.vars=c("Ind", "Pop"), variable.name = "Ancestry", value.name="Fraction")
head (mdat3)
mdat2 <- melt (Q2, id.vars=c("Ind", "Pop"), variable.name = "Ancestry", value.name="Fraction")
head (mdat2)
coloresClust <- c("royalblue3", "black", "gray80")
p3 <- ggplot(mdat3, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p3
p2 <- ggplot(mdat2, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p2
library(cowplot)
plot_grid(p2, p3, labels=c("A", "B"), ncol=1, nrow=2)
####graficamos valores de k
K <-read.table(file = "Kerror_trn.txt", header = TRUE)
####graficamos valores de k
K1 <-read.table(file = "Pleco_Reducido_Kerror_trn.txt", header = TRUE)
## grafico k
k10 <- K1[2,]
k10
k10 <- k10[-2,]
####we plot values of k
K1 <-read.table(file = "Pleco_Reducido_Kerror_trn.txt", header = TRUE)
## grafico k
k10 <- K1[2,]
k10
k1 <- k1[-2,]
K1 <- K1[-2,]
K1
K1 <- rbind(K1, k10)
plot(K1, xlab ="K value", ylim = "CV error")
plot(K1$V4, xlab ="K value", ylim = "CV error")
plot(K1$V4, xlab ="K value", ylab = "CV error")
plot(a$V4, xlab ="K value", ylab = "CV error")
plot(K1$V4, xlab ="K value", ylab = "CV error")
plot(K1$X0.71351, xlab ="K value", ylim = "CV error")
## grafico k
k10 <- K[2,]
k10
K <- K[-2,]
K
K <- rbind(K, k10)
plot(K, xlab ="K value", ylim = "CV error")
plot(K$valor, xlab ="K value", ylim = "CV error")
plot(K$valor, xlab ="K value", ylab = "CV error")
plot(K$valor, xlab ="K value", ylab = "CV error")
plot(K, pch= 19, xlab= "K value", ylab= "CV error")
plot(K$valor, xlab ="K value", ylab = "CV error")
ggplot(K, aes(x=K, y=valor))+
  geom_line(colour= "black") +
  geom_point(size=2, shape= 21, fill="white", colour= "red") +
  theme_minimal()
plot(K, pch= 19, xlab= "K value", ylab= "CV error")
pdf ("Barplo_K2_Pop_Fig3.pdf", width = 7, height = 3.5)
p3 <- ggplot(mdat3, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p3
p2 <- ggplot(mdat2, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p2
pdf ("Barplo_K2_Pop_Fig3.pdf", width = 7, height = 3.5)
p2
dev.off()
pdf ("Barplo_K3_Pop_Fig3.pdf", width = 7, height = 3.5)
p3
dev.off()
pdf ("Barplo_K2_K3_Pop_Fig3.pdf", width = 7, height = 3.5)
plot_grid(p2, p3, labels=c("A", "B"), ncol=1, nrow=2)
dev.off()
######### lo que sigue de lo anteriror antes de esto#####
colnames(Q2) <- c("Cluster 1", "Cluster 2")
Q2
colnames(Q3) <- c("Cluster 1", "Cluster 2", "Cluster 3")

DatPob <-read.table(file = "popmap1", header = FALSE)
DatPob
colnames(DatPob) <- c("Ind", "Pop")
colnames(DatPob)

Q2 <- data.frame(DatPob$Ind, DatPob$Pop, Q2)
head (Q2)
colnames(Q2) <- c("Ind", "Pop", "Cluster 1", "Cluster 2")

Q3 <- data.frame(DatPob$Ind, DatPob$Pop, Q3)
head (Q3)
colnames(Q3) <- c("Ind", "Pop", "Cluster 1", "Cluster 2", "Cluster 3")

mdat3 <- melt (Q3, id.vars=c("Ind", "Pop"), variable.name = "Ancestry", value.name="Fraction")
head (mdat3)

mdat2 <- melt (Q2, id.vars=c("Ind", "Pop"), variable.name = "Ancestry", value.name="Fraction")
head (mdat2)

coloresClust <- c("royalblue3", "black", "gray80")

p3 <- ggplot(mdat3, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p3

p2 <- ggplot(mdat2, aes(x=Ind, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ Pop, drop=TRUE, space="free", scales="free") +
  scale_fill_manual(values = coloresClust) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.text = element_text(hjust=0.5, size=8)) +
  theme(axis.title = element_blank())
p2

library(cowplot)
plot_grid(p2, p3, labels=c("A", "B"), ncol=1, nrow=2)

