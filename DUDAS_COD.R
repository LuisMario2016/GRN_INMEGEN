install.packages("renv")
library(renv)
setwd("C:\\Users\\luis_\\OneDrive\\Documentos\\INMEGEN")
renv::init()

renv::install("devtools")
library(devtools)

## instalando Summariized 
devtools::install_github("Bioconductor/SummarizedExperiment", ref = "RELEASE_3_14")

## instalando TCGA
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

## instalando ven diagram
renv::install("VennDiagram@1.6.20")

library(SummarizedExperiment)#1.22.0
library(TCGAbiolinks)#2.20.1
library(VennDiagram)#1.6.20

mthyltn <-  GDCquery(project = "TCGA-STAD",
                     data.category = "DNA Methylation",
                     platform="Illumina Human Methylation 450")
mthyltn=getResults(mthyltn)
i=substr(mthyltn$cases,1,19)

xprssn <- GDCquery(project = "TCGA-STAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type="STAR - Counts")
xprssn=getResults(xprssn)
j=substr(xprssn$cases,1,19)

mirnas <- GDCquery(project = "TCGA-STAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification")
mirnas=getResults(mirnas)
k=substr(mirnas$cases,1,19)

##############CONCOURRENT MEASURES########################
sapply(list(i,j,k),function(x) length(unique(x)))

#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)

samples=data.frame(cbind(samples,
                         sapply(samples,function(x) 
                           unique(as.character(mthyltn$sample_type[i==x])))))
colnames(samples)[2]="tissue"
samples$patient=substr(samples$samples,1,12)

subtypes=TCGAquery_subtype(tumor="stad")#subtype per patient
sum(samples$patient%in%subtypes$patient)
samples=samples[samples$patient%in%subtypes$patient,]

table(subtypes$Lauren.Class[subtypes$patient%in%samples$patient])

samples=samples[which(!samples$patient%in%subtypes$patient[
  subtypes$Lauren.Class=="NA"]|
    samples$tissue=="Solid Tissue Normal"),]

#subtype per sample
samples$subtype=sapply(samples$patient,function(x) 
  subtypes$Lauren.Class[subtypes$patient==x])
samples$subtype[samples$tissue=="Solid Tissue Normal"]="Normal"
table(samples$subtype)

write.table(samples,"subtype.tsv",sep='\t',quote=F,row.names=F)

#####################################################################
x<-head(subtype$samples,2)
xprssn <- GDCquery(project = "TCGA-STAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type="STAR - Counts", barcode = x)

GDCdownload(xprssn)
expre=GDCprepare(xprssn,summarizedExperiment=F)


mthyltn <-  GDCquery(project = "TCGA-STAD",
                     data.category = "DNA Methylation",
                     platform="Illumina Human Methylation 450",
                     data.type = "Methylation Beta Value",
                     barcode=subtype$samples)
GDCdownload(mthyltn)
