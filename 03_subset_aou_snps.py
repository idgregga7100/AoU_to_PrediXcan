#AoU Subsetting SNPs notebook 02_get_snp_intersection
#tl;dr: final path ANYONE ELSE NEEDS is gs://fc-secure-e659d40c-61b8-4bef-9728-5c017d2a3453/data/plink/predixcan_acaf_threshold.chr*.p*

#copied from dr w's script
import os
import subprocess
import numpy as np
import pandas as pd
name_of_file_in_bucket = 'predixcan_models_varids-effallele.txt'

# get the bucket name
my_bucket = os.getenv('WORKSPACE_BUCKET')

# copy csv file from the bucket to the current working space
os.system(f"gsutil cp '{my_bucket}/data/{name_of_file_in_bucket}' .")

print(f'[INFO] {name_of_file_in_bucket} is successfully downloaded into your working space')

#file has chr:pos:a1:a2,eff_allele. will need to potentially flip a1:a2
#supposedly aou varid format is chr:pos:a2:a1?
f=open('predixcan_models_varids-effallele.txt','r').readlines()
for i in (0,len(f)-1):
    f[i]=f[i].rstrip('\n') #apparently this did squat

#first build a dictionary to store chrposa1a2 with effect allele
#keys: varID, values: eff_allele
effallele = dict()
for line in f:
    row=line.split(',')
    if row[0]=='varID':
        continue
    if row[0]=='':
        continue
    row[0]=row[0].rstrip(':b38')
    effallele[row[0]]=row[1]
  
#need to just search by chr:pos and then check alleles
chrpos={}
for snp in effallele.keys():
    id=snp.split(':')
    a=id[0:2]
    b=id[2:4]
    chrpos[':'.join(a)]=':'.join(b)

#print to check
print(list(effallele.keys())[0:1])
print(list(effallele.values())[:1])
print(list(chrpos.keys())[0:1])
print(list(chrpos.values())[0:1])

#impossible to tell strand flips
badstrands=['A:T','T:A','C:G','G:C']

#make new .bim files with rsids and make list of rsids to extract from bed/bim/fam
!mkdir plink
plinkdir = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/"
for i in range(1,23):
    i = str(i)
    outlist = open("plink/chr" + i + "_keeplist","w") #list of SNPs to keep (for later --extract)
    #cp .bim to workspace
    os.system(f"gsutil -u $GOOGLE_PROJECT -m cp '{plinkdir}acaf_threshold.chr{i}.bim' plink/.")
    #count how many added to keeplist
    x=0
    y=0
    z=0
    with open("plink/acaf_threshold.chr" + i + ".bim") as f:
        for line in f:
            arr = line.strip().split()
            c = arr[0]
            p = arr[3]
            snpid = arr[1]
            id=snpid.split(':')
            a=id[0:2]
            b=id[2:4]
            ida=':'.join(a)
            idb=':'.join(b)
            
            #get snps in the effallele dict
            if ida in chrpos:
                if idb in badstrands:
                    z=z+1
                    continue
                #print(chrpos[ida],idb)
                if chrpos[ida]==idb:
                    outlist.write(snpid + '\n')
                    x=x+1
                else:
                    #print(chrpos[ida])
                    alleles=chrpos[ida].split(':')
                    flipped=alleles[::-1]
                    flipped=':'.join(flipped)
                    #print('flipped',flipped,idb)
                    if flipped==idb:
                        outlist.write(snpid+'\n')
                        y=y+1
                    else:
                        #print(snpid,chrpos[ida])
                        continue
            else:
                continue
    #check how many snps got added to the list, if at least 90% then we're good, if not then we need allele flip
    print('matched',x,'flipped-matched',y,'ambiguous skipped',z)
    #okay yeah need allele flip, around 1-2% for chr20-22
    #JUST KIDDING i mean the flip is fine to include but IT'S 2% OF THE WHOLE GENOME ON CHR 22 oh my god
    #so what i would actually have to do is check by only the chr22 snps in my list. hh. it's fine
    #aaaand i can't at this moment figure out how to do that. so. in R 69040 chr22 snps
    outlist.close()

#cp bed/bim/fam, extract needed snps and rm bed/bim/fam
#will need to increase space to run all chr's
plinkdir = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/"
for i in range(1,23):
    i = str(i)
    #cp .fam to workspace
    os.system(f"gsutil -u $GOOGLE_PROJECT -m cp '{plinkdir}acaf_threshold.chr{i}.fam' plink/.")
    #cp .bed to workspace (SLOW!)
    os.system(f"gsutil -u $GOOGLE_PROJECT -m cp '{plinkdir}acaf_threshold.chr{i}.bed' plink/.")
    #run plink2 command to only keep snps in models in pgen format
    os.system(f"plink2 --bfile plink/acaf_threshold.chr{i} --extract plink/chr{i}_keeplist --make-pgen --out plink/predixcan_acaf_threshold.chr{i}")
    #remove and list duplicte rsid SNPs
    os.system(f"plink2 --pfile plink/predixcan_acaf_threshold.chr{i} --rm-dup exclude-mismatch list --make-pgen --out plink/chr{i}")
    #replace pgen/pvar/psam files with those with duplicate rsids removed
    os.system(f"mv plink/chr{i}.pgen plink/predixcan_acaf_threshold.chr{i}.pgen")
    os.system(f"mv plink/chr{i}.pvar plink/predixcan_acaf_threshold.chr{i}.pvar")
    os.system(f"mv plink/chr{i}.psam plink/predixcan_acaf_threshold.chr{i}.psam")   
    #rm large bed/bim/fam files
    os.system(f"rm plink/acaf_threshold.chr{i}.*")
    #cp plink2 output to bucket so others can use (-o parallelizes for files > 150MB to speed up)
    os.system(f"gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp plink/predixcan_acaf_threshold.chr{i}* {my_bucket}/data/plink/")

