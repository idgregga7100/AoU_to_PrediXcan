#all unique snps across the lists from each set of models
#in aou: chr#:pos:a2:a1 essentially a bim file of all the snps is the goal
#going to have to be careful with ref/alt alleles to make sure we have the beta right

library(data.table)
library(dplyr)
library(stringr)
"%&%" = function(a,b) paste(a,b,sep="")

#concatenated models
files<-list.files('/home/isabelle/model_snps',pattern='snps.csv')
key<-data.frame(matrix(ncol=2))
colnames(key)<-c('varID','eff_allele')
for(file in files){
  f<-fread('/home/isabelle/model_snps/'%&%file)
  list<-select(f,varID,eff_allele)
  list<-unique(list)
  varid<-pull(list,varID)
  if(str_detect(varid[1],'_')){
    varid<-str_replace_all(varid,'_',':')
  }
  final<-cbind(as.data.frame(varid),list$eff_allele)
  colnames(final)<-c('varID','eff_allele')
  key<-rbind(key,final)
}

#so all the varids, and the effect allele so we can confirm
key<-unique(key)

fwrite(key,'predixcan_models_varids-effallele.txt',quote=F,row.names=F)
