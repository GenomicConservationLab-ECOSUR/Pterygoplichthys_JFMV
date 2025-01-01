######### First we performed Pearson's correlation ###########
### Libraries used
library(tidyverse)
library(PerformanceAnalytics)
library(apaTables)
library(psych)
library(corrr)
library(corrplot)
library(palmerpenguins)

read.csv("data.csv",header = TRUE,sep = ",")-> fq
attach(fq)
View(fq)
str(fq)

corrplot.mixed(fq)


######## Principal Component Analysis###########
pca <- princomp(fq) ##### generate PCA
summary(pca, loadings = T, cutoff = 0.3)

library(factoextra) #### see the eigenvalues of each component
eig_val <- get_eigenvalue(pca)
eig_val

fviz_eig(pca, addlabels = T, ylim = c(0, 50)) ### view eigenvalues in graphic mode
fviz_eig(pca, choice = c("eigenvalue"), addlabels = T, ylim = c(0, 3))

fviz_pca_biplot(pca, repel = F, col.var = "black", col.ind = "gray") ####generate a biplot

fviz_contrib(pca, choice = "var", axes = 1, top = 10)
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
fviz_contrib(pca, choice = "var", axes = 3, top = 10)

######### ANOVA ##################################################
library(tidyverse)
library(ggpubr)
library(rstatix)
library(PMCMRplus)
#### we plot to see the distribution of the data
##before plotting the data, convert the data from chr to factor

fq$sample-site<-factor(fq$sample-site)
str(fq)
##### we plot
plot(x = fq$sample-site, y = fq$Temperature,
     xlab = "Sample site", ylab = "Temperature °C", 
     col = c("gray", "gray", "gray"))

##### we calculate basic statistics ####################

fq %>%
  group_by(sample-site) %>%
  get_summary_stats(Temperature, type = "mean_sd")

model <- lm(Temperature ~ sample-site, data = fq) ##### we generated a linear model
ggqqplot(residuals(model)) ### generate a graph of the residuals of the model
shapiro_test(residuals(model)) ### we test the normality of residues

##### we can calculate normality by group
fq %>%
  group_by(sample-site) %>%
  shapiro_test(Temperature)

####### we calculate homogeneity of variances

bartlett.test(fq$Temperature ~ fq$sample-site) ###using bartlett

fq %>% levene_test(Temperature ~ sample-site) #### using Levene when any of the data is not normal

res.aov <-fq %>% anova_test(Tr ~ Sitio) #### we run the anova
res.aov

###### if there are significant differences we run poshoc

pwc <- fq %>% tukey_hsd(Temperature ~ sample-site)
pwc ### to see the results

###### we plot the results of the anova

ggboxplot(fq, x = "sample-site", y = "Temperature") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))



