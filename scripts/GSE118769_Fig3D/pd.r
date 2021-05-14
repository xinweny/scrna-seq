#Computing Pearson's correlation coefficient between cells pairwise distance genelist
#Data is already log-normalsied TPM values
#setup5.Robj is a processed single cell normalised matrix from Seuart
#setup5_genelist.txt is a list genes belonging of Robust geneset 

load("setup5.Robj")
library(Seurat)
library(dplyr)
library(Matrix)
library(Hmisc)
library(pheatmap)
library(gdata)
library(vioplot)
# dm_pos and dm_neg is based upon 820 genelist
#genepos<-read.table("dm_pos.txt",header=TRUE,quote="")
#geneneg<-read.table("dm_neg.txt",header=TRUE,quote="")
#geneall<-read.table("dm_all.txt",header=TRUE,quote="")
#genecommon<-read.table("common_all.txt",header=TRUE,quote="")
genelist<-read.table("setup5_genelist.txt",header=TRUE,quote="")
pos<-as.matrix(crist@data[paste(genelist[,1]),1:3525])
neg<-as.matrix(crist@data[paste(genelist[,1]),3526:7360])

pos_corr<-rcorr(pos,type="spearman")
neg_corr<-rcorr(neg,type="spearman")


# Finding elements of matrix and computing cell-to-cell distance
# d = sqrt(1-u/2)

upper_pos<-upperTriangle(pos_corr$r)
upper_neg<-upperTriangle(neg_corr$r)

d_pos<-sqrt((1-upper_pos)/2)
d_neg<-sqrt((1-upper_neg)/2)

names<-c("Kat2a-WT","Kat2a-NULL")
vioplot(d_neg,d_pos,names=names,col="orange")
title(main="Transcriptional heterogeneity between genotypes",ylab="Pairwise distance measure")

