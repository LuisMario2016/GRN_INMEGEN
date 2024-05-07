#!/usr/bin/env Rscript
######BEFORE concatena.R coz
#var distributions should be similar if the planned 
#approach relies on correlation networks [Tarazona2020]
#Rscript "/Users/biomario/Documents/INMEGEN/mfarnormi.R" "Diffuse" "Intestinal" "Mixed" "Unknown"

library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)

omic = commandArgs(trailingOnly=TRUE)
omic = "Diffuse" ## change acording your cancer types one by one

data=fread(paste(omic,"mtrx",sep='.'))
data=t(as.matrix(data[,2:ncol(data)],rownames=data$V1))
#factominer eats individuals in rows & variables in columns

####THE NORMALIZATION
pcs=PCA(data, scale.unit = TRUE, ncp = 2, graph = F)
#print("eigenvalue > 1 indicates the PC accounts for more variance than original variables")
#print(sum(pcs$eig[,1]>1))
w=pcs$svd$vs[1]
normi=t(data)/w
write.table(normi,paste(omic,"eigenNormi",sep='.'),sep='\t',quote=F)

