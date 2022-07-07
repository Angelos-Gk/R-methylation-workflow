#=============================LIBRARIES====================================================
#load the necessary libraries
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
#write.csv(dataframe,"name.csv")
#=============================READ AND SET THE DATA========================================

#set datadir
dataDirectory <-"C:/Users/angelosgk2/Desktop/analyses/Data/idat-2926"
#check the contents of the directory
list.files(dataDirectory) 
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
# read in the raw data from the IDAT files
targets <- read.metharray.sheet(dataDirectory, pattern="E-MTAB-2926-SampleSheet-trimmed.csv")
targets
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet
# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

#=============================SAMPLE QUALITY TESTING=======================================

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,1))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8,ylim=c(0,0.1), ylab="",names.arg=targets$Sample_Name)
abline(h=0.03,col="red")
title(ylab="Mean detection p-values", line=3.3, cex.lab=1.2)
legend(x= 1,y=0.1,"top", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
#The minfi qcReport function generates many other useful quality control plots.
#Generally, samples that look poor based on mean detection p-value will also look 
#poor using other metrics and it is usually advisable to exclude them from 
#further analysis.
# qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group,#use dev.off() in console if warning pops up
#         pdf="qcReport.pdf")
# remove poor quality samples
keep <- colMeans(detP) < 0.02
rgSet <- rgSet[,keep]
rgSet 
# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)
#=============================NORMALIZATION=================================================

#DEFINITION : removing any source of variation that is not related to biology but rather to 
#technical limitations, such as dye bias or batch effect

#For Infinium HumanMethylation450, within-array normalization concerns three main points: 
#background correction, color bias (or dye bias) adjustment and Infinium I/II-type bias correction
# background and dye corrections are done with Funnorm/Quantile functions!!!!!
#!!!!!!note:  probes on Chromosome X or Y should be filtered out to eliminate the impact 
# of sex on differential methylated analysis

mSetSqF <-preprocessFunnorm(rgSet) 
mSetSqQ <- preprocessQuantile(rgSet)
mSetSq <- preprocessIllumina(rgSet, bg.correct = TRUE, normalize = "controls",reference = 1)
mSetSqN <- preprocessNoob(rgSet, offset = 15, dyeCorr = TRUE, verbose = FALSE, dyeMethod="reference")
mSetSqS <- preprocessSWAN(rgSet, mSet = NULL, verbose = FALSE)
#bgSet <- bgcorrect.illumina(rgSet)
#conSet <- normalize.illumina.control(rgSet, reference = 1)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
#The density plots show the distribution of the beta values for each sample before 
#and after normalisation.
# visualise what the data looks like before and after normalisation
par(mfrow=c(2,3))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
#legend("top", legend = levels(factor(targets$Sample_Group)),
#       text.col=brewer.pal(5,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized (illumina method)", legend=FALSE)
densityPlot(getBeta(mSetSqF), sampGroups=targets$Sample_Group,
            main="Normalized (funnorm method)", legend=FALSE)
densityPlot(getBeta(mSetSqQ), sampGroups=targets$Sample_Group,
            main="Normalized (quantile method)", legend=FALSE)
densityPlot(getBeta(mSetSqN), sampGroups=targets$Sample_Group,
            main="Normalized (noob  method)", legend=FALSE)
densityPlot(getBeta(mSetSqS), sampGroups=targets$Sample_Group,
            main="Normalized (SWAN method)", legend=FALSE)
#legend("top", legend = levels(factor(targets$Sample_Group)),
#       text.col=brewer.pal(5,"Dark2"))

#delete the other norm objects to free-up RAM
rm(mSetSq,mSetSqF,mSetSqQ,mSetSqS)
#=============================DATA EXPLORATION=============================================

# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSqN), top=1000, gene.selection="common",pch= 16
        col=pal[factor(targets$Sample_Group)])
legend(x=-2,y=2.5, inset=c(0,0), xpd=TRUE, legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.75)

plotMDS(getM(mSetSqN), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Name)])
legend("topright", legend=levels(factor(targets$Sample_Name)), text.col=pal,
       bg="white", cex=0.65)
# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSqN), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend(x=-6,y=-0.5, legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.8, bg="white")

plotMDS(getM(mSetSqN), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend(x=-1,y=-0.5, legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.8, bg="white")

plotMDS(getM(mSetSqN), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend(x=-2,y=-0.75, legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.8, bg="white")

#=============================FILTERING====================================================
#preprocessIllumina/SWAN/Noob return a MethylSet while preprocessQuantile/Funnorm return a Genomic
#RatioSet, the dropLociWithSnps function cant work with MethylSet so we first have to use the 
#mapToGenome function in order to "convert" the MethylSet
a<-mapToGenome(mSetSqN)
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(a),rownames(detP)),] 
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(a) 
table(keep)
mSetSqFlt <- a[keep,]
mSetSqFlt
# sexfilter (if your data includes males and females, remove probes on the sex chromosomes)
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in%
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
# snpfilter (remove probes with SNPs at CpG site)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
# xhybfilter (exclude cross reactive probes) 
xReactiveProbes <- read.csv(file=paste('C:/Users/angelosgk2/Documents/R/R-4.0.2/library/methylationArrayAnalysis/extdata',
                                       '48639-non-specific-probes-Illumina450k.csv',
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt


#=============================DATA EXPLORATION==============================================

# plot the previous mds plots to see the changes
par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend(x=-2,y=2.5, legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)])
legend("topright", legend=levels(factor(targets$Sample_Name)), text.col=pal,
       cex=0.7, bg="white")
# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(1,3))
legend("bottomright", legend=levels(factor(targets$Sample_Name)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(2,3))
legend(x=-1,y=-0.5, legend=levels(factor(targets$Sample_Name)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(3,4))
legend(x=-2,y=1.5, legend=levels(factor(targets$Sample_Name)), text.col=pal,
       cex=0.7, bg="white")
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:15])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:15])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend(x=0.2,y=4, legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"),cex=0.6)
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"),cex=0.5)

#===============PROBE-WISE DIFERRENTIAL METHYLATION ANALYSIS=================================

# this is the factor of interest
cellType <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
# use the above to create a design matrix
design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- c(levels(cellType))
design
# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(ABC_DLBCL-gastric_DLBCL,
                            ABC_DLBCL-GC_DLBCL,
                            ABC_DLBCL-healthy_control,
                            gastric_DLBCL-GC_DLBCL,
                            gastric_DLBCL-healthy_control,
                            GC_DLBCL-healthy_control,
                            levels=design)
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.01
summary(decideTests(fit2,p.value=0.005))
#write.csv([enter name of dataframe here],file = file.choose(new = T))
# get the table of results for the first contrast 
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs1 <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs1)
DMPs2 <- topTable(fit2, num=Inf, coef=2, genelist=ann450kSub)
head(DMPs2)
DMPs3 <- topTable(fit2, num=Inf, coef=3, genelist=ann450kSub)
head(DMPs3)
DMPs4 <- topTable(fit2, num=Inf, coef=4, genelist=ann450kSub)
head(DMPs4)
DMPs5 <- topTable(fit2, num=Inf, coef=5, genelist=ann450kSub)
head(DMPs5)
DMPs6 <- topTable(fit2, num=Inf, coef=6, genelist=ann450kSub)
head(DMPs6)

#write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)
# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(1,2))
sapply(rownames(DMPs1)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})
par(mfrow=c(1,2))
sapply(rownames(DMPs2)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})
par(mfrow=c(1,2))
sapply(rownames(DMPs3)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})
par(mfrow=c(1,2))
sapply(rownames(DMPs4)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})
par(mfrow=c(1,2))
sapply(rownames(DMPs5)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})
par(mfrow=c(1,2))
sapply(rownames(DMPs6)[1:2], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")})

#=====================DIFFERENTIAL METHYLATION ANALYSIS OF REGIONS======================

myAnnotation1 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "ABC_DLBCL - gastric_DLBCL", arraytype = "450K",fdr = 0.01)
str(myAnnotation1)

myAnnotation2 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "ABC_DLBCL - GC_DLBCL", arraytype = "450K",fdr = 0.01)

myAnnotation3 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "ABC_DLBCL - healthy_control", arraytype = "450K",fdr = 0.01)

myAnnotation4 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "gastric_DLBCL - GC_DLBCL", arraytype = "450K",fdr = 0.01)

myAnnotation5 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "gastric_DLBCL - healthy_control", arraytype = "450K",fdr = 0.01)

myAnnotation6 <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "GC_DLBCL - healthy_control", arraytype = "450K",fdr = 0.01)

DMRs1 <- dmrcate(myAnnotation1, lambda=1000, C=2,min.cpgs = 10)
results.ranges1 <- extractRanges(DMRs1)
#how many regions?
length(results.ranges1)
#info about the top 6 regions
head(results.ranges1)

DMRs2 <- dmrcate(myAnnotation2, lambda=1000, C=2,min.cpgs = 10)
results.ranges2 <- extractRanges(DMRs2)
#how many regions?
length(results.ranges2)
#info about the top 6 regions
head(results.ranges2)

DMRs3 <- dmrcate(myAnnotation3, lambda=1000, C=2,min.cpgs = 10)
results.ranges3 <- extractRanges(DMRs3)
#how many regions?
length(results.ranges3)
#info about the top 6 regions
head(results.ranges3)

DMRs4 <- dmrcate(myAnnotation4, lambda=1000, C=2,min.cpgs = 10)
results.ranges4 <- extractRanges(DMRs4)
#how many regions?
length(results.ranges4)
#info about the top 6 regions
head(results.ranges4)

DMRs5 <- dmrcate(myAnnotation5, lambda=1000, C=2,min.cpgs = 10)
results.ranges5 <- extractRanges(DMRs5)
#how many regions?
length(results.ranges5)
#info about the top 6 regions
head(results.ranges5)

DMRs6 <- dmrcate(myAnnotation6, lambda=1000, C=2,min.cpgs = 10)
results.ranges6 <- extractRanges(DMRs6)
#how many regions?
length(results.ranges6)
#info about the top 6 regions
head(results.ranges6)
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
# draw the plot for the top DMR of each contrast

par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges1, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")


par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges2, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")


par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges3, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges4, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges5, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges6, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

#=====================GENE ONTOLOGY TESTING==============================================
# Get the significant CpG sites for each contrast at less than 1%FDR 
sigCpGs1 <- DMPs1$Name[DMPs1$adj.P.Val<0.01]
sigCpGs2 <- DMPs2$Name[DMPs2$adj.P.Val<0.01]
sigCpGs3 <- DMPs3$Name[DMPs3$adj.P.Val<0.01]
sigCpGs4 <- DMPs4$Name[DMPs4$adj.P.Val<0.01]
sigCpGs5 <- DMPs5$Name[DMPs5$adj.P.Val<0.01]
sigCpGs6 <- DMPs6$Name[DMPs6$adj.P.Val<0.01]
# First 10 significant CpGs of each contrast
sigCpGs1[1:10]
sigCpGs2[1:10]
sigCpGs3[1:10]
sigCpGs4[1:10]
sigCpGs5[1:10]
sigCpGs6[1:10]
# Total number of significant CpGs of each contrast at 1% FDR
length(sigCpGs1)
length(sigCpGs2)
length(sigCpGs3)
length(sigCpGs4)
length(sigCpGs5)
length(sigCpGs6)
# Get all the CpG sites used in the analysis to form the background
all1 <- DMPs1$Name
all2 <- DMPs2$Name
all3 <- DMPs3$Name
all4 <- DMPs4$Name
all5 <- DMPs5$Name
all6 <- DMPs6$Name
# Total number of CpG sites tested
length(all1)
length(all2)
length(all3)
length(all4)
length(all5)
length(all6)
#gometh (GO)
par(mfrow=c(1,1))
gst1 <- gometh(sig.cpg=sigCpGs1, all.cpg=all1, plot.bias=TRUE)
par(mfrow=c(1,1))
gst2 <- gometh(sig.cpg=sigCpGs2, all.cpg=all2, plot.bias=TRUE)
par(mfrow=c(1,1))
gst3 <- gometh(sig.cpg=sigCpGs3, all.cpg=all3, plot.bias=TRUE)
par(mfrow=c(1,1))
gst4 <- gometh(sig.cpg=sigCpGs4, all.cpg=all4, plot.bias=TRUE)
par(mfrow=c(1,1))
gst5 <- gometh(sig.cpg=sigCpGs5, all.cpg=all5, plot.bias=TRUE)
par(mfrow=c(1,1))
gst6 <- gometh(sig.cpg=sigCpGs6, all.cpg=all6, plot.bias=TRUE)
# Top 10 GO categories
topGSA(gst1, number=100)
topGSA(gst2, number=10)
topGSA(gst3, number=10)
topGSA(gst4, number=10)
topGSA(gst5, number=10)
topGSA(gst6, number=10)
#gometh (KEGG)
par(mfrow=c(1,1))
gst1k <- gometh(collection = "KEGG",sig.cpg=sigCpGs1, all.cpg=all1, plot.bias=TRUE)
par(mfrow=c(1,1))
gst2k <- gometh(collection = "KEGG",sig.cpg=sigCpGs2, all.cpg=all2, plot.bias=TRUE)
par(mfrow=c(1,1))
gst3k <- gometh(collection = "KEGG",sig.cpg=sigCpGs3, all.cpg=all3, plot.bias=TRUE)
par(mfrow=c(1,1))
gst4k <- gometh(collection = "KEGG",sig.cpg=sigCpGs4, all.cpg=all4, plot.bias=TRUE)
par(mfrow=c(1,1))
gst5k <- gometh(collection = "KEGG",sig.cpg=sigCpGs5, all.cpg=all5, plot.bias=TRUE)
par(mfrow=c(1,1))
gst6k <- gometh(collection = "KEGG",sig.cpg=sigCpGs6, all.cpg=all6, plot.bias=TRUE)
# Top 10 KEGG categories
topGSA(gst1k, number=100)
topGSA(gst2k, number=10)
topGSA(gst3k, number=10)
topGSA(gst4k, number=10)
topGSA(gst5k, number=10)
topGSA(gst6k, number=10)
# load Broad human curated (C2) gene sets
load(paste("C:/Users/angelosgk2/Desktop/R things/rscripts/minfi/2926 trimmed/human_c2_v5p2.rdata",sep="/"))
# perform the gene set test(s)
gsa1 <- gsameth(sig.cpg=sigCpGs1, all.cpg=all1, collection=Hs.c2,sig.genes = TRUE)
gsa2 <- gsameth(sig.cpg=sigCpGs2, all.cpg=all2, collection=Hs.c2,sig.genes = TRUE)
gsa3 <- gsameth(sig.cpg=sigCpGs3, all.cpg=all3, collection=Hs.c2,sig.genes = TRUE)
gsa4 <- gsameth(sig.cpg=sigCpGs4, all.cpg=all4, collection=Hs.c2,sig.genes = TRUE)
gsa5 <- gsameth(sig.cpg=sigCpGs5, all.cpg=all5, collection=Hs.c2,sig.genes = TRUE)
gsa6 <- gsameth(sig.cpg=sigCpGs6, all.cpg=all6, collection=Hs.c2,sig.genes = TRUE)




# top 10 gene sets
topGSA(gsa1, number=10)
topGSA(gsa2, number=10)
topGSA(gsa3, number=10)
topGSA(gsa4, number=10)
topGSA(gsa5, number=10)
topGSA(gsa6, number=10)



sessionInfo()
