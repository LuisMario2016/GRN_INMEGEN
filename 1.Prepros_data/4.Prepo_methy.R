setwd("/Users/biomario/Documents/INMEGEN/Prepo_methy")
library(TCGAbiolinks)
library(data.table)

subtypeSTAD=read.table("subtype.tsv",header=T,sep='\t')
subtypeSTAD$subtype[is.na(subtypeSTAD$subtype)] <- "Unknown"
#get the data

mthyltn <-  GDCquery(project = "TCGA-STAD",
                     data.category = "DNA Methylation",
                     data.type = "Methylation Beta Value",
                     platform="Illumina Human Methylation 450",
                     barcode=subtypeSTAD$samples)

#GDCdownload(mthyltn)
files=list.files("GDCdata",full.names=T,recursive=T)

###############################################
#BiocManager::install("sesameData", force=T)
#BiocManager::install("sesame", force=T)
library(sesameData)
library(sesame)
mthyltn2=GDCprepare(mthyltn)
methy <- mthyltn2@assays@data[[1]]
#dim(methy)
#[1] 485577    75
#write.table(methy,"methy.tsv",sep='\t',quote=F)

#
# methy_2=BiocGenerics::do.call(cbind,pbapply::pbsapply(files,function(x) 
#   fread(x,select=2,stringsAsFactors=F)))
#######drop probes with too many na########################################
#before noise-prone filter or the process will be slooow

#subtype to duplicate
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtypeSTAD[c(which(!subtypeSTAD$samples%in%j),
                          as.numeric(sapply(which(subtypeSTAD$samples%in%j),rep,2))),1:7]
#needed coz names are not equal to expression data but barcodes do
designMethy$barcode=unlist(sapply(designMethy$samples,function(x)
  colnames(methy)[i==x][1]))
designMethy$barcode[designMethy$samples%in%j]=rev(colnames(methy)[
  which(i%in%j)])
designMethy=designMethy[order(match(designMethy$barcode,
                                    colnames(methy))),]
total=table(subtypeSTAD$subtype)
#NA per subtype
nas=lapply(names(total),function(x) 
  rowSums(is.na(methy[,designMethy$subtype==x])))
#keep probes with NA in less than 25% of samples of all subtypes
i=unique(unlist(lapply(1:2,function(x) which(nas[[x]]<total[x]*.25))))
methy=methy[i,]
#dim(methy)
#[1] 417322    75
############filter noise-prone probes##############################
#get probes annotation
annot=fread(files[1])
colnames(annot)[1]="probe"
annot=annot[annot$probe%in%rownames(methy),]
annot=annot[order(match(annot$probe,rownames(methy))),]
sum(annot$probe==rownames(methy)) #order is the same

#sex chrs should be dropped if mixed gender
#since this is not the case, keep chrX
#https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
#methy=methy[!annot$Chromosome%in%c("chrX","chrY"),]
#dim(methy)

#drop probes with ambigous mapping
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
methy_=methy[annot$Chromosome!="*",]
nrow(methy_)
# [1] 0

#polymorphisms may affect DNAm measurements
methy_grep=methy[grep("rs",rownames(methy),invert=T),]
nrow(methy_grep)
#[1] 417322 all values for elements that do not match.

#######impute missing data########################################
library(doParallel)
#BiocManager::install("pbapply", force = T)
library(impute)
library(pbapply)

#separate per subtype
methy=lapply(names(total),function(x) 
  methy[,designMethy$subtype==x])
nas=lapply(methy,function(x) rowSums(is.na(x)))
sapply(nas,max)
# [1] 70  3
sapply(nas,function(x) sum(x>0))#cell to estimate
# [1] 44762  1298

#impute missing values
no_cores=detectCores()-1
registerDoParallel(cores=no_cores)
cl=makeCluster(no_cores,type="FORK")
#takes a looong time (Mac M1 8gb 8 cores, takes 2:30 hours)
methy=parLapply(cl,methy,function(x) 
  impute.knn(x,k=15,maxp=nrow(x))$data)#maxp forces impute, even if too many NAs
#Buuren and Groothuis suggested that 15 variables are generally 
#sufficient for imputation
stopCluster(cl)
#alternative
#methy=pbapply::pblapply(methy,function(x) impute.knn(x,k=15,maxp=nrow(x))$data)

write.table(do.call(cbind,methy),"methyNormi.tsv",sep='\t',quote=F)
methydf <- read.table("methyNormi.tsv",header=T,sep='\t')

#######beta to M values########################################
#Beta values have severe heteroscedasticity for highly methylated or
#unmethylated CpG sites. M-values provide much better performance in
#terms of detection rate and true positive rate for both highly 
#methylated and unmethylated CpG sites
library(ggplot2)

#check beta distributions
temp=lapply(methy,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:2,function(x) 
  cbind(temp[[x]],levels(as.factor(designMethy$subtype))[x]))))
colnames(temp)=c("beta","subtype")
temp$beta=as.numeric(as.character(temp$beta))
png("betaDistr.png")
ggplot(temp,aes(x=beta))+
  geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
               alpha=0.15)+theme(text=element_text(size=18))
dev.off()

#transform beta to M as in:
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
mval=function(beta){log2(beta/(1-beta))}
m=pbapply::pblapply(methy,function(x) apply(x,c(1,2),mval))
#check M distributions
temp=lapply(m,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:2,function(x) 
  cbind(temp[[x]],levels(as.factor(designMethy$subtype))[x]))))
colnames(temp)=c("M","subtype")
temp$M=as.numeric(as.character(temp$M))
png("MDistr.png")
ggplot(temp,aes(x=M))+
  geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
               alpha=0.1)+theme(text=element_text(size=18))
dev.off()
# Warning message: The dot-dot notation (`..scaled..`) was deprecated in ggplot2 3.4.0.
# ℹ Please use `after_stat(scaled)` instead.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 
#######average duplicates########################################
m=do.call(cbind,m)
i=substr(colnames(m),1,19)
j=i[duplicated(i)]
prefi=m[,!i%in%j]
duplis=m[,i%in%j]

temp=do.call(cbind,lapply(j,function(x) 
  rowMeans(duplis[,substr(colnames(duplis),1,19)%in%x])))
#identify samples with barcode 

###temp es de una sola dimnsión por lo que no podemos aplicarle la funcion###
#colnames(temp)=j
################
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
write.table(final,"methyM.tsv",sep='\t',quote=F)
#2511M when tar.gz

###########check final data#################################
# BiocManager::install("NOISeq")
library(NOISeq)


#u're dragging batch effects, use them as covariates when DM
noiseqData = NOISeq::readData(data = final, factor=subtypeSTAD)
myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = T, logtransf = T)
#logtransf=F fails & norm=F flattens points
png("PCA_methy_global.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "subtype")})
dev.off()
png("PCA_methy_global_stage.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "tissue")})
dev.off()
### check if we need race and gender? i dont think so

# png("PCA_methy_global_race.png")
# print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
#                   factor = "race")})
# dev.off()