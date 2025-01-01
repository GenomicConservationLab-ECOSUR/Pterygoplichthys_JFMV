############ run wilcoxon for samples
library(tidyverse)
library(rstatix)
library(ggpubr)

read.csv("p1.csv",h=T,sep = ",")-> p1
###### resumen estadistico
p1 %>%
  group_by(group) %>%
  get_summary_stats(He, type = "median_iqr")

###### vizualizacion
bxp <- ggboxplot(
  p1, x = "group", y = "He", 
  ylab = "Expected heterozygosity (He)", xlab = "Calling method", add = "jitter"
)
bxp

##### calculation

stat.test <- p1 %>% 
  wilcox_test(He ~ group) %>%
  add_significance()
stat.test

##### effect size
p1 %>% wilcox_effsize(He ~ group)

##### graphical report
stat.test <- stat.test %>% add_xy_position(x = "group")
bxp + 
  stat_pvalue_manual(stat.test, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test, detailed = TRUE))
