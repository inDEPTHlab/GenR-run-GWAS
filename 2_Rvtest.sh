#!/bin/bash

# Generation R analyses
# Within GENA, navigate to your project folder and run this script by copying it in the command window
# For information on rvtest, see http://zhanxw.github.io/rvtests/
# SCRIPT NO 2, what this script does:
#    - Main analysis: GWAS:  continuous standardised phenotype ~ SNP + sex + child age + first 10 genetic principal components
#    - Secondary analysis 1: continuous standardised phenotype ~ SNP + child age + first 10 genetic principal components AMONG GIRLS ONLY
#    - Secondary analysis 2: continuous standardised phenotype ~ SNP + child age + first 10 genetic principal components AMONG BOYS ONLY

################################################################# MAIN ANALYSIS 1Y ##########################################################################
for c in {1..22}
do
    if [ -f ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz ]
    then
        echo "Chromosome $c"
        qsub-bin-long rvtest --inVcf ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz --dosage DS --pheno sourcefiles/phenotype.ped --pheno-name phenotype_z --covar sourcefiles/covariates.ped --covar-name sex,age,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --out results/mainanalysis/chr${c} --single score --noweb
    fi
done

################ SECONDARY 1: GIRLS ONLY ################
for c in {1..22}
do
    if [ -f ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz ]
    then
        echo "Chromosome $c (girls)"
        qsub-bin-long rvtest --inVcf ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz --dosage DS --pheno sourcefiles/phenotype_girls.ped --pheno-name phenotype_z --covar sourcefiles/covariates.ped --covar-name age,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --out results/girlsonly/chr${c} --single score --noweb
    fi
done

################ SECONDARY 1: BOYS ONLY #################
for c in {1..22}
do
    if [ -f ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz ]
    then
        echo "Chromosome $c (boys)"
        qsub-bin-long rvtest --inVcf ~/GENR3/Imputed/1000G_PhaseIIIv5/chr${c}.dose.vcf.gz --dosage DS --pheno sourcefiles/phenotype_boys.ped --pheno-name phenotype_z --covar sourcefiles/covariates.ped --covar-name age,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --out results/boysonly/chr${c} --single score --noweb
    fi
done

################ X-chromosome ################

qsub-bin-long rvtest --inVcf ~/GENR3/Imputed/1000G_PhaseIIIv5/chrX.no.auto_male.dose.vcf.gz --dosage DS --pheno sourcefiles/phenotype.ped --pheno-name phenotype_z --covar sourcefiles/covariates.ped --covar-name age,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --out results/mainanalysis/chrXmale --single score --noweb

qsub-bin-long rvtest --inVcf ~/GENR3/Imputed/1000G_PhaseIIIv5/chrX.no.auto_female.dose.vcf.gz --dosage DS --pheno sourcefiles/phenotype.ped --pheno-name phenotype_z --covar sourcefiles/covariates.ped --covar-name age,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --out results/mainanalysis/chrXfemale --single score --noweb
