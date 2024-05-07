setwd("/Users/biomario/Documents/INMEGEN/Prepo_miRNAS")

library(TCGAbiolinks)
library(biomaRt)

subtypeSTAD=read.table("subtype.tsv",header=T,sep='\t')
subtypeSTAD[subtypeSTAD == "Unknown"] <- NA
#subtypeSTAD$subtype[is.na(subtypeSTAD$subtype)] <- "Unknown"
#designExp$subtype[is.na(designExp$subtype)] <- "Unknown"
#get the data
mirnas <- GDCquery(project = "TCGA-STAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification",
                   barcode=subtypeSTAD$samples)

#Genome of reference: hg38
#https://api.gdc.cancer.gov/data/683367de-81c9-408c-85fd-2391f3e537ee

#GDCdownload(mirnas)
mir=GDCprepare(mirnas)
rownames(mir)=mir$miRNA_ID
mir=mir[,grep("read_count",colnames(mir))]
colnames(mir)=gsub("read_count_","",colnames(mir))
dim(mir)
dim(mir)
# 1881  242
write.table(mir,"miRNAseq.tsv",sep='\t',quote=F)

#subtype to duplicates
i=substr(colnames(mir),1,19)
j=i[duplicated(i)]
designExp=subtypeSTAD[c(which(!subtypeSTAD$samples%in%j),
                    as.numeric(sapply(which(subtypeSTAD$samples%in%j),rep,2))),]
designExp=designExp[order(match(designExp$samples,substr(colnames(mir),1,19))),]
designExp$barcode=colnames(mir)

mart=useMart("ensembl",host="https://oct2022.archive.ensembl.org", 
             dataset = "hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "mirbase_id",
                             "start_position","end_position"),
              mart=mart)
myannot=myannot[myannot$mirbase_id%in%rownames(mir),]
#there should not be a length bias
myannot$length=abs(myannot$end_position-myannot$start_position)
#discard duplicated entries with the same %CpG
sum(!duplicated(myannot[,2:3]))
# [1] 1863
myannot=myannot[!duplicated(myannot[,2:3]),]
#there're duplicates with slightly different %CpG & position
# temp=myannot[myannot$mirbase_id%in%myannot$mirbase_id[duplicated(myannot$mirbase_id)],]
# summary(colSums(sapply(unique(temp$mirbase_id),function(x)
# temp$percentage_gene_gc_content[temp$mirbase_id==x])*c(1,-1)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.4600 -0.4200  0.3300  0.1647  0.5150  1.7400 
#choose 1 of the duplicates randomly
myannotAlt=myannot[duplicated(myannot$mirbase_id),]
myannot=myannot[!duplicated(myannot$mirbase_id),]
#if the GC bias aint fixed you CAN NOT compare among miRNAs

##################CHECK BIASES########################################################
library(NOISeq)

noiseqData = NOISeq::readData(data = mir, factor=designExp,
                              gc=myannot[,c(3,2)],length=myannot[,c(3,6)])
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio",factor = "subtype")#check low counts
#distributions
png("miROri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:2)
dev.off()
#counts
png("miRcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:2)
dev.off()
png("miRlowCountThres.png")
hist(rowMeans(edgeR::cpm(mir,log=T)),ylab="miRNA",
     xlab="mean of log CPM",col="gray",xlim=c(-5,20))
dev.off()

#check length & GC bias
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias",
                   factor = "subtype")
png("miRGCbiasOri.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias",
                 factor = "subtype")
png("lengthbiasOri.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(mylenBias, samples = x))
dev.off()

myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = FALSE, 
                    logtransf = FALSE)#check batches
png("miRPCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
           factor = "subtype")
dev.off()
mycd = dat(noiseqData, type = "cd", norm = FALSE)#check if normalizations is needed
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Warning: 366 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: TCGA-BR-8592-01A-11R-2402-13"
# FAILED 
# 241
png("miRcdOri.png")
explo.plot(mycd,samples=1:10)
dev.off()
#################SOLVE BIASES######################################################
#filter low counts
FilteredMatrix = filtered.data(mir, factor = "subtype",
                               norm = FALSE, method = 1, cpm = 0)
#Filtering out low count features...
#227 features are to be kept for differential expression analysis with filtering method 1
# of  the  samples 
#temp=lapply(unique(designExp$subtype),function(x)   #this was commented
#  mir[,colnames(mir)%in%designExp$barcode[designExp$subtype==x]])
#temp1=names(which(table(unlist(sapply(temp,function(x) rownames(x)[rowSums(x>=5)>=ncol(x)*.25])))==5))
#length(temp)
# [1] 4

#FilteredMatrixAlt=mir[rowSums(mir)>0,]

#TMM, UQ, median & DESEq are similar [10.1186/gb-2010-11-3-r25]
#TMM, UQ are the best [10.1093/bib/bbv019] 
myTMM=tmm(FilteredMatrix,lc=0)
myTMMAlt=tmm(FilteredMatrixAlt,lc=0)
noiseqData = NOISeq::readData(data = myTMM, factors=designExp)
#noiseqData = readData(data = myTMMAlt, factors=designExp)
mycdTMM = NOISeq::dat(noiseqData, type = "cd", norm = T)
#mycdTMMAlt = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#1    240 

myUQ=uqua(FilteredMatrix,lc=0)
noiseqData = NOISeq::readData(data = myUQ, factors=designExp)
mycdUQ = NOISeq::dat(noiseqData, type = "cd", norm = T)
table(mycdUQ@dat$DiagnosticTest[,  "Diagnostic Test"])
# FAILED PASSED 
# 61    180  

library(EDASeq)
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(FilteredMatrix),
  phenoData=data.frame(designExp,row.names=designExp$barcode))

norm.counts <- betweenLaneNormalization(mydataEDA,
                                        which = "median", offset = FALSE)
noiseqData = NOISeq::readData(data = assayData(norm.counts)$normalizedCounts,
                              factors=designExp)
mycdMedian = NOISeq::dat(noiseqData, type = "cd", norm = T)
table(mycdMedian@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
# 7    234

library("DESeq2")
#deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix,
#                                                 conditions=designExp))
#myDESEQ=counts(deseqFactors,normalized=T)
#noiseqData = NOISeq::readData(data = myDESEQ, factors=designExp)
#mycdDESEQ = NOISeq::dat(noiseqData, type = "cd", norm = T)
#table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   707    101 

#png("miRcd_final.png")
#explo.plot(mycdMedian,samples=1:10)
#dev.off()
rownames(designExp)<-designExp$barcode
dds <- DESeqDataSetFromMatrix(countData = FilteredMatrix, colData = designExp, design = ~ 1)

# estimar factores de escala de tamaño
dds <- estimateSizeFactors(dds)

write.table(FilteredMatrix, "FiltermatrixSTAD.tsv", sep ="\t", quote=F)

myDESEQ=counts(dds,normalized=T)
noiseqData = NOISeq::readData(data = myDESEQ, factors=designExp)
mycdDESEQ = NOISeq::dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
# FAILED PASSED 
# 6      235 

png("miRcd_final.png")
explo.plot(mycdMedian,samples=1:10)
dev.off()

#############################SOLVE BATCH EFFECT#######################################################
noiseqData = NOISeq::readData(data = assayData(norm.counts)$normalizedCounts,
                              factors=designExp)
myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = T,logtransf=F)#log=F or points'll fall in angle
png("miRPCA_preArsyn.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "subtype")
dev.off()
nobatch=ARSyNseq(noiseqData, factor = "subtype", batch = F,
                 norm = "n",  logtransf=F)#log=F or dots'll collapse

#############################FINAL QUALITY CHECK#######################################################
noiseqData = NOISeq::readData(data = exprs(nobatch),
                              factors=designExp)
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio",factor = "subtype")#check low counts
png("miRFinal.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:2)
dev.off()
png("miRcountsFinal.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:2)
dev.off()
myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = T,logtransf=F)#check batches
png("miRPCA_Final.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "subtype")
dev.off()
#############################RESOLVE DUPLICATES & SAVE##################################################
miRfinal=assayData(norm.counts)$normalizedCounts
#get duplicated index
i=designExp$samples[duplicated(designExp$samples)]
#get barcodes per sample
i=lapply(i,function(x) designExp$barcode[designExp$samples==x])
#same duplicates than in RNAseq, difference's still the plate
#separate duplicates
duplis=miRfinal[,colnames(miRfinal)%in%unlist(i)]
prefi=miRfinal[,!colnames(miRfinal)%in%unlist(i)]
#average duplicates
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
colnames(prefi)=substr(colnames(prefi),1,19)
#checar y revisar duplis y temp
#colnames(duplis)
#filas_seleccionadas <- c("TCGA-FP-A4BF-01A-12R-A360-13", "TCGA-FP-A4BF-01A-12R-A24L-13", "TCGA-HU-8249-01A-11R-2343-13",
#                         +                          "TCGA-HU-8249-01A-11R-A360-13", "TCGA-BR-A44U-01A-11R-A360-13", "TCGA-BR-A44U-01A-11R-A24L-13",
#                         +                          "TCGA-CD-A48A-01A-12R-A360-13", "TCGA-CD-A48A-01A-12R-A24L-13", "TCGA-HU-A4GQ-01A-11R-A360-13",
#                         +                          "TCGA-HU-A4GQ-01A-11R-A252-13", "TCGA-HU-A4GD-01A-11R-A252-13", "TCGA-HU-A4GD-01A-11R-A360-13")
#informacion_seleccionada <- subset(designExp, rownames(designExp) %in% filas_seleccionadas)
#subtypeSTAD$samples %in% informacion_seleccionada$samples
#summary(subtypeSTAD$samples %in% informacion_seleccionada$samples)
#colnames(temp)<-subtype_seleccionado$samples


#joint matrices
final=cbind(prefi,temp)
dim(final)
# [1] 244  75
final=final[,order(match(colnames(final),subtypeSTAD$barcode))]
#subtypeSTAD$samples %in% colnames(final)
write.table(final,"miRNAseqNormi.tsv",sep='\t',quote=F)
