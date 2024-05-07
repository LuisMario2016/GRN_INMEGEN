
setwd("/Users/biomario/Documents/INMEGEN/")
library(data.table)
subtypeSTAD=read.table("subtype.tsv",header=T,sep='\t')
subtypeSTAD$subtype[is.na(subtypeSTAD$subtype)] <- "Unknown"
# revert comments for unnormalized data
# expre=fread("RNAseqnormalized.tsv")
# miR=fread("miRNAseqNormi.tsv")
# methy=fread("methyM.tsv")
files=list.files("./multiomics",full.names=T)
#files=files[grep("eigenNormi",files)]
data=lapply(files,fread)
#expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
#miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
#methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)

data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(data)=gsub("./multiomics/","",files)
names(data)=gsub(".tsv","",names(data))
names(data)=gsub("RNAseqnormalized","transcripts",
                 gsub("miRNAseqNormi","miRNAs",gsub("methyM","CpGs",names(data))))
print(sapply(data,dim))
#print(sapply(data,function(x) head(rownames(x))))

#choose methy order
#subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#expre=expre[,order(match(colnames(expre),subtype$samples))]
#miR=miR[,order(match(colnames(miR),subtype$samples))]
subtype=subtypeSTAD[order(match(subtypeSTAD$samples,colnames(data$CpGs))),]
data[2:3]=lapply(data[2:3],function(x)
  x[,order(match(colnames(x),subtype$samples))])
print(names(data))

subtypeSTAD$subtype <- as.factor(subtypeSTAD$subtype)

#data per subtype
concatenated=lapply(levels(subtypeSTAD$subtype),function(x) 
  list(CpGs=data$CpGs[,subtypeSTAD$subtype==x],
       transcripts=data$transcripts[,subtypeSTAD$subtype==x],
       miRNA=data$miRNAs[,subtypeSTAD$subtype==x]))

names(concatenated)=levels(subtypeSTAD$subtype)


##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
print(sapply(concatenated,dim))
#Diffuse Intestinal Mixed Unknown
#[1,]  418087     418087 418087  418087
#[2,]      58        152     19       7
lapply(1:4,function(x) write.table(concatenated[[x]],
                                   paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))
#adjust 1:4 acording your cancer types

################################ for 1.6
#run in terminal
# Rscript 1_6mfanormi.R LUSC