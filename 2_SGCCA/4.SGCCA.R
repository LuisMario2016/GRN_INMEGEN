#!/usr/bin/env Rscript
setwd("/Users/biomario/Documents/INMEGEN/")
########################PARAMETERS & PACKAGES
args=commandArgs(trailingOnly=TRUE)
subtype="Mixed"


library(igraph)
library(mixOmics)
library(data.table)
########################DATA
data=fread(paste("/Users/biomario/Documents/INMEGEN/",
                 paste(subtype,"eigenNormi",sep='.'),sep='/'))
data=as.matrix(data[,2:ncol(data)],rownames=data$V1)
#separate omics
first_letter <- substr(rownames(data),1 ,1)
newData = list()
for (letter in unique(first_letter)) {
  newData[[letter]] <- data[first_letter == letter,] 
}
names(newData)=c("CpGs","transcripts","miRNAs")
data <- newData
rm(newData) #remove to free memory
penalty=c(CpGs=0.24,transcripts=0.44,miRNAs=0.25)#output of ObtainPenalty.R

transposeData <- function(df){
  data = list()
  data[[1]] = t(df[[1]])
  data[[2]] = t(df[[2]])
  data[[3]] = t(df[[3]])
  names(data)=c("CpGs","transcripts","miRNAs")
  return(data)
}
########################THE SGCCA
ncomp=ncol(data$miRNAs)-1#the last comp has all loadings>0
final=wrapper.sgcca(X=transposeData(data),penalty=penalty,scale=F,
                    scheme="centroid",ncomp=ncomp, max.iter=100000)#ncomp to explain 50% of transcripts matrix according to mfa.R
#get selected features
selected=lapply(final$loadings,function(y)
  apply(y,2,function(x) x[x!=0]))
transcripts_list <- lapply(1:ncol(selected$transcripts), function(i) {
    selected$transcripts[, i]})
View(transcripts_list)
selected$transcripts <- transcripts_list
View(selected)
names(selected$transcripts) <- paste0("comp", 1:length(selected$transcripts))
selected=as.data.frame(do.call(rbind,lapply(selected,function(y)
  do.call(rbind,lapply(1:length(y),function(x)
    cbind(names(y)[x],y[[x]],names(y[[x]])))))))
colnames(selected)=c("component","final","variable")

#####PLOT LOADINGS
library(ggplot2)
library(gridExtra)

selected$omic=substr(selected$variable,1,1)
selected$omic=gsub("E","transcripts",
                   gsub("h","miRNAs",gsub("c","CpGs",selected$omic)))
selected$final=as.numeric(as.character(selected$final))

png(paste(subtype,"loadings.png",sep='-'))
ggplot(selected,aes(x=omic,y=final))+
  geom_boxplot()+ylab("loading")+theme(text=element_text(size=18))
dev.off()

###############
initial=wrapper.sgcca(X=transposeData(data),penalty=rep(1,3),scale=F,
                      scheme="centroid",ncomp=ncomp, max.iter=100000)
rbind(rowSums(do.call(rbind,initial$AVE$AVE_X)),
      rowSums(do.call(rbind,final$AVE$AVE_X)))

selected=lapply(unique(selected$omic),function(x) 
  selected[selected$omic==x,])
temp=lapply(1:3,function(y) apply(selected[[y]],1,function(x) 
  initial$loadings[[y]][x[3],x[1]]))
selected=do.call(rbind,selected)
selected$initial=unlist(temp)
plots=lapply(unique(selected$omic),function(x)
  ggplot(selected[selected$omic==x,],
         aes(y=final,x=initial))+geom_point()+ggtitle(x)+
    theme(text=element_text(size=18)))
png(paste(subtype,"loadings_change.png",sep='-'))
grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

setwd("/Users/biomario/Documents/INMEGEN/SCCCA")
write.table(selected,paste(subtype,"selected",sep='.'),sep='\t',
            quote=F,row.names=F)
