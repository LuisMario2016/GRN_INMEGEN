
setwd("/Users/biomario/Documents/INMEGEN/")
temp=read.table("penalty_search.tsv",sep='\t',header=T)

#######################################DIAGNOSTIC PLOTS
library(ggplot2)#3.3.5
library(gridExtra)#2.3

# omics as factor
temp$omic=factor(temp$omic,levels=c("CpGs","transcripts","miRNAs"))
temp=temp[order(temp$penalty),]
# remove NA's
temp <- temp[!is.na(temp$omic),]

# parse variables to numeric
temp$AVE = as.numeric(temp$AVE)
temp$nfeatures = as.numeric(temp$nfeatures)
temp$penalty = as.numeric(temp$penalty)
str(temp)

# plot AVE and number of features
png("AVE.png",width=800)
ggplot(temp,aes(y=AVE,x=as.character(round(penalty,2)),
                group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
  theme(text=element_text(size=18),
        axis.text.x = element_text(angle = 45))
dev.off()
png("nfeatures.png",width=800)
ggplot(temp,aes(y=nfeatures,x=as.character(round(penalty,2)),
                group=penalty))+geom_boxplot()+facet_wrap(~omic)+xlab("penalty")+
  theme(text=element_text(size=18),
        axis.text.x = element_text(angle = 45))+
  scale_y_continuous(trans="log10")
dev.off()
#plot median AVE vs meadian nfeatures
omics=levels(temp$omic)
omics=lapply(omics,function(x) temp[temp$omic==x,])
# omics=lapply(omics,function(x) as.data.frame(apply(x,2,as.numeric))) #this not
names(omics)=levels(temp$omic)
omics=lapply(omics,function(y) sapply(unique(y$penalty),function(x) 
  apply(y[y$penalty==x,],2,median,na.rm=T)))#better than mean?
omics=lapply(omics,function(x) as.data.frame(t(x)))
str(omics)

# parse variable to numeric (this could be a function)
omics$CpGs$AVE = as.numeric(omics$CpGs$AVE)
omics$CpGs$nfeatures = as.numeric(omics$CpGs$nfeatures)
omics$CpGs$penalty = as.numeric(omics$CpGs$penalty)

omics$transcripts$AVE = as.numeric(omics$transcripts$AVE)
omics$transcripts$nfeatures = as.numeric(omics$transcripts$nfeatures)
omics$transcripts$penalty = as.numeric(omics$transcripts$penalty)

omics$miRNAs$AVE = as.numeric(omics$miRNAs$AVE)
omics$miRNAs$nfeatures = as.numeric(omics$miRNAs$nfeatures)
omics$miRNAs$penalty = as.numeric(omics$miRNAs$penalty)

#indi plots or CpGs will determine axis
plots=lapply(1:3,function(x) ggplot(omics[[x]],
                                    aes(x=nfeatures,y=AVE,col=penalty))+geom_point()+
               ggtitle(names(omics)[x])+theme(text=element_text(size=18))+
               scale_x_continuous(trans="log10")+geom_line())
png("sparsity_search.png")
grid.arrange(plots[[1]],plots[[2]],plots[[3]])
dev.off()

# choose penalty
grid=unique(temp$penalty)
#slopes=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
  #y$AVE[x+1]-y$AVE[x]))
slopes1=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
  abs(y$AVE[x+1]-y$AVE[x])/abs(y$nfeatures[x+1]-y$nfeatures[x])))


slopes1$miRNAs[abs(slopes1$miRNAs)=="Inf"]=NA
slopes1$transcripts[abs(slopes1$transcripts)=="Inf"]=NA
penalty=sapply(slopes,function(x) grid[which.max(x)+1])
penalty1=sapply(slopes1,function(x) grid[which.max(x)+1])

# > penalty
# CpGs transcripts      miRNAs 
# 0.03        0.60        0.60 

# > penalty1
# CpGs transcripts      miRNAs 
# 0.01        0.02        0.11 


#CHOOSE penalty1

#######################################PENALTIES SUGGESTED BY PLOTS by SoL
grid=unique(temp$penalty)
#slopes=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
#	y$AVE[x+1]-y$AVE[x]))
slopes1=lapply(omics,function(y) sapply(1:(length(grid)-1),function(x) 
  abs(y$AVE[x+1]-y$AVE[x])/abs(y$nfeatures[x+1]-y$nfeatures[x])))

slopes1$miRNAs[abs(slopes1$miRNAs)=="Inf"]=NA
