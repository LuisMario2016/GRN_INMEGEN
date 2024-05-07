###########################################################
setwd("/Users/biomario/Documents/INMEGEN/")
penalty_cpgs <- 0.24
penalty_transcripts <- 0.44
penalty_mir <- 0.25

# ' Read files in correct form
# '@param dataTable String name of the file. You add .csv, .tsv, etc.
# '@return data_table data_frame with the correct form
readMtrx <- function(dataTable){
  df <- data.table::fread(dataTable)
  return(as.matrix(df[,2:ncol(df)], rownames=df$V1))
}

#correct subtype
##subtypeSTAD$subtype[is.na(subtypeSTAD$subtype)] <- "Unknown"
##write.table(subtypeSTAD, "subtype.tsv", sep = '\t', quote = F)
subtype <- readMtrx("/Users/biomario/Documents/INMEGEN/subtype.tsv")
subtype_subtype <- as.factor(subtype[,"subtype"])
summary(subtype_subtype)
library(mixOmics)
#take model descriptors<-----------------recycled
describe=function(data,pc,pt,pm){
  #subsample observations
  #data=lapply(data,function(x) x[sample(1:n,subn),])
  resus=wrapper.sgcca(data,penalty=c(pc,pt,pm),scale=T,
                      scheme="centroid")
  #get results description
  description=as.data.frame(do.call(rbind,resus$AVE$AVE_X))
  description$nfeatures=sapply(resus$loadings,function(x) sum(x!=0))
  description$omic=rownames(description)
  description$penalty=resus$penalty
  colnames(description)[1]="AVE"
  return(description)
}

transposeData <- function(df){
  data = list()
  data[[1]] = t(df[[1]])
  data[[2]] = t(df[[2]])
  data[[3]] = t(df[[3]])
  names(data)=c("CpGs","transcripts","miRNAs") # we suppose we work with this omic
  return(data)
}

library(data.table)
library(parallel)

files=list.files()
files=files[grep("eigenNormi",files)]
sizes=c(Diffuse=58,Intestinal=152,Mixed=19,Unknown=7)


###
# there is a bias-variance trade-off associated with the choice of k-fold
# cross-validation. Typically, given these considerations, one performs k-fold 
# cross-validation using k=5 or k=10, as these values have been shown empirically
# to yield test error rate estimates that suffer neither from excessively high bias
# nor from very high variance. Hastie, Tibshirani et. all 2021
## 
# we use CV with k=5 for each tumor subtype 
k = 5 # k folds
# we have problem due to we have only 5 samples from normal tissue so in each 
# iteration we choose 3 samples from normal tissue randomly taken

# load data
data=lapply(1:length(files),function(x) fread(files[x]))
names(data) = gsub(".eigenNormi", "", files)

data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))

set.seed(2) #we will to use a seed to replicate results
index = list()
for (j in 1:length(files)) {index[[j]] = list()}

## tumor subtypes indexs without replace
 #for(i in 1:length(files)){ # from 1 or 2??
   #temp_vec <- cut(sample(1:dim(data[[i]])[2]), breaks = 5, labels = F, ordered_result = T)
  # index[[i]] <- split(1:dim(data[[i]])[2], temp_vec)
 #}

# get indexs of samples for each fold
for (i in 1:length(files)) {
  len_vec <- sizes[i]
  samps <- replicate(k, sample(1:len_vec, size = round(len_vec*0.2)))
  
  if(class(samps)[1] == "integer"){ 
    # this is in case 20% of samples are less or equal to 1,in this situation
    # we get a vector numeric instead of matrix array
    for(j in 1:k){index[[i]][[j]] <- samps[j]}
  } else {
    for (j in 1:k) {
      index[[i]][[j]] <- samps[,j]
    }
  }
  
}

# get samples for each fold and add to dataframe by omics
descr <- data.frame()
for (fold in 1:k) {
  
  selected_index = list()
  for (j in 1:length(files)) {
    selected_index[[j]]=index[[j]][[fold]]
  }
  
  selectData = lapply(1:length(files), 
                      function(x) data[[x]][,selected_index[[x]]])
  selectData=do.call(cbind,selectData)
  
  #separate omics
  newData = list()
  first_letter <- substr(rownames(selectData),1 ,1)
  # > unique(first_letter)
  # [1] "c" "E" "h"
  for (letter in unique(first_letter)) {
    newData[[letter]] <- selectData[first_letter == letter,] 
  }
  names(newData)=c("CpGs","transcripts","miRNAs") #based on unique(first_letter)
  # confirm separates omics well
  #head(rownames(newData$CpGs)); tail(rownames(newData$CpGs))
  #head(rownames(newData$transcripts)); tail(rownames(newData$transcripts))
  #head(rownames(newData$miRNAs)); tail(rownames(newData$miRNAs))
  
  descr <- rbind(descr, describe(transposeData(newData),
                                 penalty_cpgs,penalty_transcripts,penalty_mir))

}

setwd("/Users/biomario/Documents/INMEGEN/")
file=paste(penalty_cpgs, penalty_transcripts, penalty_mir, "tsv", sep = '.')
write.table(descr, file, sep='\t', quote=F, row.names = F)
###########################