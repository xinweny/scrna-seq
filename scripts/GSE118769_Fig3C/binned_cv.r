#Calculating squared coefficient of variation sCV and DM
# POS 3525 NEG 3835 Total 7360 Genes 2588
# Only normalised values are used

load("setup5.Robj")
#load("setup6.Robj")
library(Seurat)
library(dplyr)
library(Matrix)
library(zoo)

genelist<-read.table("setup5_genelist.txt",header=TRUE,quote="")

crepos <- grep("1", rownames(crist@data.info), value = T)
creneg <- grep("2", rownames(crist@data.info), value = T)

CV_matrix_gcn5<-matrix(,nrow=length(genelist[,1]),ncol=9)

for ( j in 1:length(genelist[,1])){
	CV_matrix_gcn5[j,1]<-paste(genelist[j,1])
	# For CREPOS cells
	print(genelist[j,1])
	#pos<-crist@raw.data[gcn5[j],1:length(crepos)]
	pos<-crist@data[genelist[j,1],1:length(crepos)]
	CV_matrix_gcn5[j,2]<-round(mean(pos),4)
	CV_matrix_gcn5[j,3]<-round(sd(pos),4)
	CV_matrix_gcn5[j,4]<-round(((sd(pos)/mean(pos))^2),4) #squared CV
	
	# For CRENEG cells
	#neg<-crist@raw.data[gcn5[j],(length(crepos)+1):7360]
	neg<-crist@data[genelist[j,1],(length(crepos)+1):7360]
	CV_matrix_gcn5[j,6]<-round(mean(neg),4)
	CV_matrix_gcn5[j,7]<-round(sd(neg),4)
	CV_matrix_gcn5[j,8]<-round(((sd(neg)/mean(neg))^2),4) #squared CV
}

# This is stored as globalCVmatrix.txt, consists of mean, sd and CV values for CREPOS and CRENEG cells for robust geneset genes
cv<-(CV_matrix_gcn5) 

colnames(cv)<-c("gene","pos_mean","pos_sd","pos_sCV","Fano_pos","neg_mean","neg_sd","neg_sCV","Fano_neg")
pos<-cbind(cv[,2],cv[,4]) # mean and sCV
neg<-cbind(cv[,6],cv[,8]) # mean and sCV

# For plotting mean values were divided into intervals in the range (-0.5 to 0.9) to generate binned CV plot.
par(font.lab=2,font.axis=2)
at=c(1:10)
neg_cut<-(cut(log10(as.numeric(neg[,1])),breaks=c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.1,0.1,0.3,0.9)))
pos_cut<-(cut(log10(as.numeric(pos[,1])),breaks=c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.1,0.1,0.3,0.9)))
boxplot(log10(as.numeric(pos[,2]))~(pos_cut),boxwex=0.25,col="orange",xaxt="n",at = at + 0.15, ylab="log10(squared CV)", xlab= "mean expression bins" )
boxplot(log10(as.numeric(neg[,2]))~(neg_cut),boxwex=0.25,col="grey", add = T,at = at - 0.15 )
legend("topright",legend=c("Kat2a-WT","Kat2a-NULL"),col=c("grey","orange"),pch=16)
dev.print(pdf,file="sCV_global.pdf")


