#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
argms=commandArgs(trailingOnly=TRUE)
subtype="Diffuse"
iteration=10
#seed for replicate analysis
set.seed(42 + as.numeric(iteration))
#ncomp=as.numeric(args[2]]))
setwd("/Users/biomario/Documents/INMEGEN/")

#library(igraph)
library(mixOmics)
library(data.table)
########################DATA
data=fread(paste(subtype,"eigenNormi",sep='.'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)

#sample data
i=ncol(data)
data=data[,sample(1:i,round(i/2))]

#separate omics
first_letter <- substr(rownames(data),1 ,1)
newData = list()
for (letter in unique(first_letter)) {
  newData[[letter]] <- data[first_letter == letter,] 
}
names(newData)=c("CpGs","transcripts","miRNAs")
data <- newData

transposeData <- function(df){
  data = list()
  data[[1]] = t(df[[1]])
  data[[2]] = t(df[[2]])
  data[[3]] = t(df[[3]])
  names(data)=c("CpGs","transcripts","miRNAs")
  return(data)
}

########################THE SGCCA
penalty=c(CpGs=0.01,transcripts=0.02,miRNAs=0.11)#output of choose_penalty.R
ncomp=ncol(data$miRNAs)-1#the last comp has all loadings>0


final=wrapper.sgcca(X=transposeData(data),penalty=penalty,scale=F,
                    scheme="centroid",ncomp=ncomp)#ncomp to explain 50% of transcripts matrix according to mfa.R
#get selected features
selected=do.call(rbind,lapply(1:3,function(y) 
  #y covers the 3 data blocks
  do.call(rbind,lapply(1:ncomp,function(x) 
    cbind("comp"=x,
          "feature"=selectVar(final,comp=x,block=y)[[1]][[1]])))))
#selectVar output is a list of length 2,1st element is a list with name & value)
setwd("/Users/biomario/Documents/INMEGEN/SCCCA")
write.table(selected,paste(subtype,iteration,"selected",sep='.'),sep='\t',
            quote=F,row.names=F)

