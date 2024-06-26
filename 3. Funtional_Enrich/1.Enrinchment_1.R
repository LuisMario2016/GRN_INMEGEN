setwd("/datos/home/mario2016/INMEGEN")
library(tidyverse)
library(biomaRt)
library(rentrez)
library(ggplot2)
#load sgcca.R output
files=list.files()
files=files[grepl("\\.selected$", files) & !grepl("\\.subtype\\.selected$", files)]
sets=lapply(files,read_tsv)
names(sets)=gsub(".selected","",files)

##################CHECK FEATURE STABILITY####################
files=list.files()
files=files[grep("subsampled",files)]#output of joinSubsamples.R
#files <- files[-c(1, 3)] Deleted replicates if its case
stability=lapply(files,read_tsv)
names(stability)=gsub(".subtype.selected.subsampled","",files)

names(stability)==names(sets)#to be sure
stability=lapply(stability,function(x) as.data.frame(cbind(
  "variable"=x$feature,
  #sum across subsamples
  "stability"=rowSums(!is.na(x[,3:ncol(x)])))))

#add stability info to sets
subtypes=names(sets)
sets=lapply(subtypes,function(x) 
  merge(sets[[x]],stability[[x]], by=))
names(sets)=subtypes#to be sure

##################FOR THE PLOT
#join all together in a data.frame
stability=do.call(rbind,lapply(1:4,function(x) 
  as.data.frame(cbind("subtype"=names(stability)[x],
                      stability[[x]]))))
#nice omic identifiers
stability$omic=gsub("E","transcript",
                    gsub("h","miRNA",
                         gsub("c","CpG",substr(stability$variable,1,1))))
#force the order that I want
stability$omic=factor(stability$omic,levels=c("CpG",
                                              "transcript",
                                              "miRNA"))
#categorical instead of numbers
stability$tokeep=c("non stable","stable")[
  as.factor(stability$stability>=1)]#arbitrary threshold
#barplot
png("stability.png")
stability%>%count(subtype,omic,tokeep)%>%
  ggplot(aes(y=as.numeric(n),x=subtype,fill=tokeep))+
  geom_bar(stat="identity",position="fill")+facet_wrap(~omic)+
  xlab("")+ylab("%")+theme(text=element_text(size=16),
                           panel.background=element_blank(),legend.position="bottom",
                           legend.title = element_blank(),
                           axis.text.x=element_text(angle=45),
                           legend.margin=margin(-40,0,0,0))+#or the legen fall far bottom
  scale_fill_manual(values=c("gray47","brown1"))
dev.off()

#from now one use just stable features
sets=lapply(sets,function(x) x%>%filter(stability>=1))
#setwd("~/lungsquamouscells/functional_enrichment")
lapply(names(sets),function(x)
  write_tsv(x=sets[[x]],file=paste(x,"stable",sep='.')))


####################FEATURE INTERSECTIONS################### 
library(UpSetR)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)  

#up-set plot per omic & subtype of selected features
omics=unique(sets$Diffuse$omic)
features=lapply(omics,function(x) lapply(sets,function(y)
  #list features selected
  y%>%filter(omic==x)%>%distinct(variable)%>%unlist))
pdf("slctd_feats_intrsctn.pdf",width=20,height=10)
lapply(features,function(x) 
  upset(fromList(x),text.scale=3,order.by="degree"))
dev.off()

#top selected features
library(dplyr)
library(purrr)

top=lapply(sets,function(x) x%>%dplyr::count(variable,omic)
%>%arrange(desc(n))%>%group_by(omic)%>%group_map(~head(.x,5)))
top=lapply(top,function(x) do.call(rbind,x))
top=as.data.frame(do.call(rbind,lapply(1:4,function(x)
  cbind(top[[x]],subtype=names(top)[x]))))
#add gene descriptions
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
                version=105)
#https://dec2021.archive.ensembl.org
myannot=getBM(attributes = c("ensembl_gene_id","hgnc_symbol",
                             "wikigene_description"),mart=mart)
colnames(myannot)[1]="variable"


#add gene lfc
#setwd("/Users/biomario/Documents/INMEGEN/")
de=read_tsv("DE.genes.tsv") #dont have DE.genes.tsv
colnames(de)[1:2]=c("subtype","variable")
#de$subtype=gsub("_Normal","",de$subtype)
topg=merge(top,myannot,by="variable")
merge(de[,c(1:3,7)],topg,by=c("variable","subtype"))
topg <- arrange(topg, subtype, desc(n))
topg<-write_tsv(topg, "top.tsv")

#repeated features
top=top%>%pivot_wider(names_from="subtype",values_from=n)
top[rowSums(is.na(top))<4,]
