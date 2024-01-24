# AoU_to_PrediXcan

The process of going from the snps in our models on wl3 to running a GWAS in AoU and then taking the sumstats back to wl3 for PrediXcan.

/home/isabelle/model_snps

AoU workspace Subsetting SNPs

tl;dr: final path to subsetted aou snps for predixcan is gs://fc-secure-e659d40c-61b8-4bef-9728-5c017d2a3453/data/plink/predixcan_acaf_threshold.chr*.p*

After running GWAS with subsetted plink files and concatenating chrs, download final sumstats file to local computer, then use scp to copy to wheelerlab3 for running Pred/SPred. NOTE: AoU snps in variant ID format, SPred uses chrpos. i split the snp col of the sumstats files to get chrpos in R, and then SPred ran with 95% of snps in the model. Pred might have a flag to use the varid in the models, but SPred does not
