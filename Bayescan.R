#### bayescan analysis
# bayescan analysis
# first: transform the vcf file to .geste with PGDSpider
# second: we run bayescan on command line
### bayescan2 ./pleco.geste -od ./ -threats 2 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10

# we read the output file in R _fst.txt
bayescan=read.table("pleco.g_fst.txt") 

# we extract the SNP list from the vcf file in the terminal
# grep -v "#" 3059snps_103ind.vcf | cut -f 3 > id-3059snps.txt

# then import the SNP-ids into R
SNPb=read.table("id-3059snps.txt",header=FALSE)

# merge the outlier names with the results of the bayescan data frame

bayescan1=cbind(SNPb, bayescan) 

#rename the columns
colnames(bayescan1)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

#type results
write.table(bayescan1, "SNP_pleco_bayesresult.txt", quote=FALSE, sep="\t", row.names=FALSE)
## change the value of column Q_VALUE
attach(bayescan1)
class(bayescan1$Q_VALUE)

bayescan1$Q_VALUE <- as.numeric(bayescan1$Q_VALUE) 
bayescan1[bayescan1$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

#Round the values
bayescan1$LOG_PO <- (round(bayescan1$LOG_PO, 4)) 
bayescan1$Q_VALUE <- (round(bayescan1$Q_VALUE, 4)) 
bayescan1$ALPHA <- (round(bayescan1$ALPHA, 4)) 
bayescan1$FST <- (round(bayescan1$FST, 6))

## we add a colum of type of selection
bayescan1$SELECTION <- ifelse(bayescan1$ALPHA>=0&bayescan1$Q_VALUE<=0.05,"diversifying",ifelse(bayescan1$ALPHA>=0&bayescan1$Q_VALUE>0.05,"neutral","balancing")) 
bayescan1$SELECTION<- factor(bayescan1$SELECTION)
levels(bayescan1$SELECTION)

## save SNP results under positive and balancing selection
positive <- bayescan1[bayescan1$SELECTION=="diversifying",] 
neutral <- bayescan1[bayescan1$SELECTION=="neutral",] 
balancing <- bayescan1[bayescan1$SELECTION=="balancing",]

#Check the number of SNP in each category
xtabs(data=bayescan1, ~SELECTION)

#write potential snp results under selection
write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F)

# transform log Q values for plotting

range(bayescan1$Q_VALUE) 
bayescan1$LOG10_Q <- -log10(bayescan1$Q_VALUE) 

####plot with ggplot
# 1) we create the titles of the axes

x_title="Log(q-value)" 
y_title="Fst" 

# we elaborate the chart
ggplot(bayescan1,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+   
  theme_classic()

