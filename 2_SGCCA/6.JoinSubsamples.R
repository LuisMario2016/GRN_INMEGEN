#setwd("/Users/biomario/Documents/INMEGEN/")
setwd("/Users/biomario/Documents/INMEGEN/SCCCA")

library(tidyverse)

subtype=commandArgs(trailingOnly=TRUE)
subtype <- "Intestinal.10"

files=list.files()
files=files[grep(subtype,files)]
subsa=lapply(files,read_tsv)
#subsa <- lapply(files, function(file) read_tsv(file, locale = locale(encoding = "latin1")))
#final=read_tsv(paste("../",subtype,".10.selected",sep=''))
final <- read_tsv("Unknown.selected")
subsa <- subsa[-c(2, 3, 4)]

subsa=lapply(subsa,function(x) x%>%dplyr::count(subsa[[1]][["feature"]]))
names(subsa[[1]])[1] <- "feature"
subsa=reduce(subsa,full_join,by="feature")
colnames(subsa)[2:101]=paste("r",1:100,sep="")
#subsa <- lapply(subsa, function(df) {
 # colnames(df) <- c("feature", paste("r", 1:100, sep = ""))
  #return(df)
#})
final=final%>%count(variable)
colnames(final)[1]="feature"
subsa=merge(final,subsa,by="feature",all.x=T)
write_tsv(subsa,paste(subtype,"subsampled",sep='.'))
#temp=apply(subsa[,2:101],1,table)

