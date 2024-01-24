#split varid from AoU sumstats to chrpos for running SPred
library(dplyr)
library(stringr)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

sumstats<-list.files('/home/isabelle/model_snps/height_sumstats',pattern='AoU')

for (file in sumstats){
stats<-fread(file)
varid<-pull(stats,SNP)
varid<-str_split_fixed(varid,':',3)
chrpos<-str_c(varid[,1],varid[,2],sep=':')%>%as.data.frame()
stats<-cbind(chrpos,stats)
colnames(stats)<-c('chrpos','SNP','A1','A2','BETA','SE')
fwrite(stats,'chrpos-'%&%file,quote=F,row.names=F,sep='\t')
}
