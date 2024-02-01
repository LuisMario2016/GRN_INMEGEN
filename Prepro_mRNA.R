#####################################################################
library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)  
setwd("/Users/biomario/Documents/INMEGEN")
subtype <- read.delim("~/Documents/INMEGEN/subtype.tsv")
xprssn <- GDCquery(project = "TCGA-STAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type="STAR - Counts", barcode = subtype$samples)

#GDCdownload(xprssn)
expre=GDCprepare(xprssn,summarizedExperiment=F)
expreT=GDCprepare(xprssn,summarizedExperiment=T)
tempSTAD = as.matrix(expre[,2:ncol(expre)])
rownames(tempSTAD) = expre$gene_id
expre = tempSTAD

stranded_firstSTAD <- expreT@assays@data@listData[["stranded_first"]]
unstrandedSTAD <- expreT@assays@data@listData[["unstranded"]]
stranded_secondSTAD <- expreT@assays@data@listData[["stranded_second"]]
tpm_unstrandSTAD <- expreT@assays@data@listData[["tpm_unstrand"]]
fpkm_unstrandSTAD <- expreT@assays@data@listData[["fpkm_unstrand"]]
fpkm_uq_unstrandSTAD <- expreT@assays@data@listData[["fpkm_uq_unstrand"]]
# "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"

variables_ <- list(stranded_firstSTAD, unstrandedSTAD, stranded_secondSTAD, tpm_unstrandSTAD, 
                   fpkm_unstrandSTAD, fpkm_uq_unstrandSTAD)

colnames(stranded_firstSTAD) <- expreT@colData@rownames
colnames(unstrandedSTAD) <- expreT@colData@rownames
colnames(stranded_secondSTAD) <- expreT@colData@rownames
colnames(tpm_unstrandSTAD) <- expreT@colData@rownames
colnames(fpkm_unstrandSTAD) <- expreT@colData@rownames
colnames(fpkm_uq_unstrandSTAD) <- expreT@colData@rownames

write.table(unstrandedSTAD, "RNAseqSTAD.tsv", sep = '\t', quote = F)

# subtype to duplicates #one only
i = substr(colnames(unstrandedSTAD), 1, 19)
j = i[duplicated(i)]
designExpSTAD=subtype[c(which(!subtype$samples%in%j),
                        as.numeric(sapply(which(subtype$samples%in%j),rep,2))),]
designExpSTAD=designExpSTAD[order(match(designExpSTAD$samples,substr(colnames(expre),1,19))),]
designExpSTAD$barcode=colnames(fpkm_unstrandSTAD)

expre[,"gene_name"][1:5]
colnames(expre)[1:5]

add_gene_info <- function(dfSTAD){
  dfSTAD <- cbind(gene_type = expre[,"gene_type"][1:60660], dfSTAD)
  dfSTAD <- cbind(gene_name = expre[,"gene_name"][1:60660], dfSTAD)
  dfSTAD <- cbind(gene_id = rownames(expre)[1:60660], dfSTAD)
  return(dfSTAD)
}
unstrandedSTAD <- add_gene_info(unstrandedSTAD)
stranded_firstSTAD <- add_gene_info(stranded_firstSTAD)
stranded_secondSTAD <- add_gene_info(stranded_secondSTAD)
tpm_unstrandSTAD <- add_gene_info(tpm_unstrandSTAD)
fpkm_unstrandASTAD <- add_gene_info(fpkm_unstrandSTAD)
fpkm_uq_unstrandSTAD <- add_gene_info(fpkm_uq_unstrandSTAD)

dim(designExpSTAD)

# keep only tenscript id not version numbers
rownames(unstrandedSTAD) <- unstrandedSTAD[,"gene_id"]
rownames(unstrandedSTAD) <- sapply(strsplit(rownames(unstrandedSTAD), ".", fixed=T),
                                   function(x) x[1])

rownames(stranded_firstSTAD) <- rownames(unstrandedSTAD)
rownames(stranded_secondSTAD) <- rownames(unstrandedSTAD)
rownames(tpm_unstrandSTAD) <- rownames(unstrandedSTAD)
rownames(fpkm_unstrandASTAD) <- rownames(unstrandedSTAD)
rownames(fpkm_uq_unstrandSTAD) <- rownames(unstrandedSTAD)

## "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"
#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "gene_biotype",
                             "start_position","end_position","hgnc_id","hgnc_symbol"),
              filters = "ensembl_gene_id", 
              values=rownames(fpkm_unstrandASTAD),mart=mart) #its valid for every variable

myannot$length=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&
                  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=unstrandedSTAD[rownames(fpkm_uq_unstrandSTAD)%in%myannot$ensembl_gene_id,]
exprots_hgnc <- exprots_hgnc[,4:ncol(exprots_hgnc)]
dim(exprots_hgnc)

##check duplicated probes
myannot2 <- myannot[unique(rownames(myannot)),]
dim(myannot2); dim(myannot)

myannot$hgnc_id[duplicated(myannot$hgnc_id)]
# [1] "HGNC:30046" "HGNC:11582" "HGNC:33853" "HGNC:4876"  
which(myannot2$hgnc_id == "HGNC:30046"); which(myannot2$hgnc_id == "HGNC:32460")
# [1] 18566 18728
# [1] 19327 19336
which(myannot2$hgnc_id == "HGNC:33853"); which(myannot2$hgnc_id == "HGNC:4876")
# [1] 12080 19340
# [1]  7581 19375

myannot2[c(18573,18736),]; myannot2[c(16090,19260),]
myannot2[c(12080,19356),]; myannot2[c(7581,19393),]

myannot3 <- myannot2[-c(18736,12080,19260,7581),]
dim(myannot2); dim(myannot3)
# [1] 19382     8
# [1] 19378     8
length(unique(rownames(stranded_firstSTAD)))
# [1] 60616

##################CHECK BIASES########################################################
library(NOISeq)
library(edgeR)

exprots_hgnc2 <- as.data.frame(exprots_hgnc[unique(rownames(exprots_hgnc)),],)
exprots_hgnc3 <- sapply(exprots_hgnc2, as.numeric)
rownames(exprots_hgnc3) <- rownames(exprots_hgnc[unique(rownames(exprots_hgnc)),])


#format data for noiseq
noiseqData = NOISeq::readData(data = exprots_hgnc3,
                              gc = myannot[,1:2],
                              biotype = myannot[,c(1,3)],factor=designExpSTAD,
                              length=myannot[,c(1,8)])
noiseqData2 = NOISeq::readData(data = exprots_hgnc3,
                               gc = myannot3[,1:2],
                               biotype = myannot3[,c(1,3)],factor=designExpSTAD,
                               length=myannot3[,c(1,8)])

#1)check expression bias per subtype
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio", factor = "subtype")
# [1] "Warning: 200 features with 0 counts in all samples are to be removed for this analysis."
# [1] "Counts per million distributions are to be computed for:"
# [1] "LUSC"   "normal"
mycountsbio2 = NOISeq::dat(noiseqData2, type = "countsbio", factor = "subtype")
# [1] "Warning: 200 features with 0 counts in all samples are to be removed for this analysis."
# [1] "Counts per million distributions are to be computed for:"
# [1] "LUSC"   "normal"

#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio2, plottype = "boxplot", samples = 1:2)
dev.off()
#2)check for low count genes
png("lowcountsOri.png")
explo.plot(mycountsbio2, plottype = "barplot", samples = 1:2)
dev.off()
png("lowCountThres.png")
hist(rowMeans(cpm(exprots_hgnc3,log=T)),ylab="genes",
     xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias
#each sample s is compared to a reference r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). 
#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
mycd = dat(noiseqData2, type = "cd", norm = FALSE) #slow
#[1] "Warning: 200 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: TCGA-18-4721-01A-01R-1443-07"
# [1] "Diagnostic test: FAILED. Normalization is required to correct this bias."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
# FAILED PASSED 
# 63     11 
png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:ncol(exprots_hgnc3),10))
dev.off()

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%) the exp,ression depends on the feature
# Asignar un valor por defecto (por ejemplo, "Unknown") a las muestras con NA
noiseqData2@phenoData@data[["subtype"]][is.na(noiseqData2@phenoData@data[["subtype"]])] <- "Unknown"

myGCcontent <- dat(noiseqData2, type = "GCbias", factor = "subtype")
png("GCbiasOri.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2, function(x) explo.plot(myGCcontent, samples = x))
dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- dat(noiseqData2, k = 0, type = "lengthbias",
                 factor = "subtype")
png("lengthbiasOri.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(mylenBias, samples = x))
dev.off()
#BUT, since the gene has the same length in all your samples, there is no need to divide by the gene length

#5) check for batch effect
myPCA = dat(noiseqData2, type = "PCA", norm = F, logtransf = F)
png("PCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "subtype")
dev.off()

#################SOLVE BIASES######################################################
library(EDASeq)

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 
countMatrixFiltered = filtered.data(exprots_hgnc3, factor = "subtype",
                                    norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#10943 features are to be kept for differential expression analysis with filtering method 1
myannot3=myannot3[myannot3$ensembl_gene_id%in%rownames(countMatrixFiltered),]
myannot=myannot[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered),]


#all names must match

# Obtén las filas que son diferentes entre los dos data frames
rows_to_remove <- setdiff(rownames(countMatrixFiltered), myannot3$ensembl_gene_id)

# Filtra countMatrixFiltered para retener solo las filas que no están en rows_to_remove
countMatrixFiltered <- countMatrixFiltered[setdiff(rownames(countMatrixFiltered), rows_to_remove), ]

mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot3,row.names=myannot3$ensembl_gene_id),
  phenoData=data.frame(designExpSTAD,row.names=designExpSTAD$barcode))


countMatrixFiltered <- countMatrixFiltered[myannot3$ensembl_gene_id, ]

#countMatrixFiltered <- countMatrixFiltered[rownames(countMatrixFiltered) %in% myannot$ensembl_gene_id, ]
# Filtrar myannot para incluir solo las filas presentes en countMatrixFiltered
myannot <- myannot[myannot$ensembl_gene_id %in% rownames(countMatrixFiltered), ]

mydataEDA2 <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(designExpSTAD,row.names=designExpSTAD$barcode))
#order for less bias
gcFull <- withinLaneNormalization(mydataEDA2, 
                                  "percentage_gene_gc_content", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 
fullfullTMM <-NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
norm.counts <- betweenLaneNormalization(normCounts(lFull),
                                        which = "median", offset = FALSE)

noiseqData = NOISeq::readData(data = fullfullTMM, factors=designExpSTAD)
#cd has to preceed ARSyN or won't work
mycd=NOISeq::dat(noiseqData,type="cd",norm=TRUE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"]) #sometimes change values
#FAILED PASSED 
#  3    232 

#############################SOLVE BATCH EFFECT#######################################################
# Asignar un valor por defecto (por ejemplo, "Unknown") a las muestras con NA
noiseqData@phenoData@data[["subtype"]][is.na(noiseqData@phenoData@data[["subtype"]])] <- "Unknown"

myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = F)
png("preArsyn.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "subtype")
dev.off()

ffTMMARSyn=ARSyNseq(noiseqData, factor = "subtype", batch = F,
                    norm = "n",  logtransf = T)

myPCA = dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = T)
png("postArsyn22.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "subtype")
dev.off()

explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "race")
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "gender")
# #############################FINAL QUALITY CHECK#######################################################
designExpSTAD$subtype[is.na(designExpSTAD$subtype)] <- "Unknown"

noiseqData <- NOISeq::readData(data = fullfullTMM, gc = myannot[,1:2],
                               biotype = myannot[,c(1,3)],factor=designExpSTAD,
                               length=myannot[,c(1,8)])
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio", factor = "subtype",
                          norm=T)

png("CountsFinal.png")
explo.plot(mycountsbio, plottype = "boxplot",samples=1:2)
dev.off()
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias",
                   factor = "subtype",norm=T)
png("GCbiasFinal.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias",
                 factor = "subtype",norm=T)
png("lengthbiasFinal.png",width=1000)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(mylenBias, samples = x))
dev.off()

#############################RESOLVE DUPLICATES & SAVE##################################################
#get duplicates
i=designExpSTAD$samples[duplicated(designExpSTAD$samples)]
#get sample barcode per sample
i=lapply(i,function(x) designExpSTAD$barcode[designExpSTAD$samples==x])
#separate duplicates
final=fullfullTMM
duplis=final[,colnames(final)%in%unlist(i)]
prefi=final[,!colnames(final)%in%unlist(i)]
#average duplicates ### THERE ARE NOT DUPLICATES EN LUSC
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
#identify samples with barcode 
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
dim(final)
# [1] 11768    75
final=final[,order(match(colnames(final),subtype$samples))]
dim(myannot3)
write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)





