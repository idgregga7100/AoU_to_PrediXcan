#run this script to concatenate all the chromosome GWAS results

#install needed R packages
lapply(c('dplyr','R.utils'),
       function(pkg) { if(! pkg %in% installed.packages()) {  install.packages(pkg) } } )
#library(qqman)
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

#read in gwas results per pop and plot
for(pop in c("afr","amr","eur")){
    res = fread("gwas/GWAS_AoU_testing_" %&% pop %&% "_chr1.scaled_ht.glm.linear.gz")
    for(i in 2:22){
    chrres = fread("gwas/GWAS_AoU_testing_" %&% pop %&% "_chr" %&% i %&% ".scaled_ht.glm.linear.gz")
    res = rbind(res,chrres)
    }
    #subset to ADD rows and plot
    res = dplyr::filter(res,TEST=="ADD")
    #remove NA's before plotting
    res = dplyr::filter(res, !is.na(P)) |> mutate(P=as.numeric(P))
    #print(manhattan(res,chr="#CHROM",bp="POS",snp="ID",main=pop,col=c("orange","black")))
    #print(qq(res$P,main=pop))
    #write sumstats to file for prs-csx
    print(dim(res))
    #select needed data for prs-csx
    sumstats = dplyr::select(res, ID, ALT, REF, BETA, SE)
    colnames(sumstats) = c("SNP","A1","A2","BETA","SE")
    print(head(sumstats))
    fwrite(sumstats, "gwas/GWAS_AoU_testing_" %&% pop %&% "_scaled_ht.sumstats.txt",row.names = FALSE,sep='\t')
}

#download final sumstats file to local computer, then use scp to copy to wheelerlab3 for running Pred/SPred
#note: AoU snps in variant ID format, SPred uses chrpos. i split the snp col of the sumstats files 
#to get chrpos in R, and then SPred ran with 95% of snps in the model. Pred might have a flag to use the varid
#in the models, but SPred does not
