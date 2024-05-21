getwd()
library(edgeR)
options(digits=3)
Coinfection.targets
install.packages("ggplot2")
library(ggplot2)
ggplot(Coinfection.rawCount) +
  geom_histogram(aes(x = Ha1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
png("./results/count distribution.png", res=300, height=1800, width=1800)
dev.off()
write.csv(Coinfection.rawCount, file="./results/Coinfection.rawCounts.csv")
infection.normCPM <- cpm(calcNormFactors(Coinfection.orig))
dim(infection.normCPM)
write.csv(infection.normCPM, file="./results/infection.normCPM.csv")
infection.filtered <- rowSums(cpm(Coinfection.orig)>1) >=3
table(infection.filtered)
Coinfection.orig$samples$lib.size
Infection <- Coinfection.orig[infection.filtered,]
colSums(Infection$counts)
Infection$samples$lib.size <- colSums(Infection$counts)
Infection$samples
Infection = calcNormFactors(Infection)
Infection$samples
Infection.filtered.normCPM <-cpm(calcNormFactors(Infection))
write.csv(Infection.filtered.normCPM, file="./results/Infection.filtered.normCPM.csv")
group<-factor(c('Ha','Ha','Ha',"Ctr","Ctr","Ctr"))

Infection.design <- model.matrix(~group)   
rownames(Infection.design)<-colnames(Infection$counts)
Infection.design
plotMDS(Infection, main="MDS plot of RNA-Seq", labels=colnames(Infection$counts))

png("./results/plotMDS.Infection.png", res=300, height=1800, width=1800)
Infection <- estimateGLMCommonDisp(Infection, Infection.design)
Infection <- estimateGLMTrendedDisp(Infection, Infection.design)
Infection <- estimateGLMTagwiseDisp(Infection, Infection.design)
plotMeanVar(Infection, show.tagwise.vars=T,NBline=T)
plotBCV(Infection)
Infection.fit <- glmFit(Infection, Infection.design)

colnames(Infection.fit)
lrt.Ha_vs_Ctr <- glmLRT(Infection.fit, coef=2)
t1<-topTags(lrt.Ha_vs_Ctr, n=nrow(Infection))
head(t1$table)
summary(decideTests(lrt.Ha_vs_Ctr, adjust.method="BH", p.value=0.05))
nrow(subset(topTags(lrt.Ha_vs_Ctr, n=586)$table,  logFC > 0))
lrt.Ha_vs_Ctr_UP <- subset(topTags(lrt.Ha_vs_Ctr, n=586)$table, logFC > 0)
nrow(subset(topTags(lrt.Ha_vs_Ctr, n=586)$table,  logFC < 0))
lrt.Ha_vs_Ctr_DW <- subset(topTags(lrt.Ha_vs_Ctr, n=586)$table, logFC < 0)
DEtags.lrt.Ha_vs_Ctr <- rownames(Infection)[as.logical(decideTests(lrt.Ha_vs_Ctr, adjust.method="BH", p.value=0.05))]
write.csv(lrt.Ha_vs_Ctr_UP, file="./results/lrt.Ha_vs_Ctr_UP.csv")

write.csv(lrt.Ha_vs_Ctr_DW, file="./results/lrt.Ha_vs_Ctr_DW.csv")
Infection.colHavsCtr = rep('grey55', nrow(Infection))
Infection.colHavsCtr[lrt.Ha_vs_Ctr$table$PValue < 0.05 & lrt.Ha_vs_Ctr$table$logFC >0 ] <- "red"
Infection.colHavsCtr[lrt.Ha_vs_Ctr$table$PValue < 0.05 & lrt.Ha_vs_Ctr$table$logFC <0 ] <- "blue"
par(omi=c(0.1,0.1,0.1,0.1), las=1, cex=0.5, mgp=c(3,1,0), cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
plotSmear(lrt.Ha_vs_Ctr, de.tags=DEtags.lrt.Ha_vs_Ctr, xlab="log-counts per million (logCPM)", ylab="log2-fold change (log2FC)", main="Ha infection compared to Control", pch=19, cex=0.4, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, ylim=c(-7,7), yaxs="i")
abline(h=c(-1,1),col="dodgerblue")
par(omi=c(0.1,0.1,0.1,0.1), las=1, cex=0.5, mgp=c(3,1,0), cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
plotSmear(lrt.Ha_vs_Ctr, xlab="log-counts per million (logCPM)", ylab="log2-fold change (log2FC)", main="a infection compared to Control", smearWidth=0.5, pch=21, cex=0.4, deCol="red", col=Infection.colHavsCtr, ylim=c(-7,7), yaxs="i")
abline(h=c(-1,1),col="dodgerblue")
png("./results/plotSmear.InfectionRNAseq.png", res=300, height=1800, width=1800)


