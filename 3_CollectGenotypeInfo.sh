#!/bin/bash

# Generation R analyses
# Within GENA, navigate to your project folder and run this script by copying it in the command window
# SCRIPT NO 3, what this script does:
#   - Get file with the Hardy-Weinburgh Equilibrium (HWE) information
#   - Get the R2 information (concerning imputation quality)
#   - Compress the GWAS results
#   - Merge the GWAS result for each chromosome into just one file

# Get HWE
plink --bfile ~/GENR3/Genotyped/GENR3-2012 --hardy --out sourcefiles/GenR # genotyped only
gzip sourcefiles/GenR.hwe

# Get R squared, info score etc 
zcat ~/GENR3/Imputed/1000G_PhaseIIIv5/chr22.info.gz | head -n1 > sourcefiles/rsq.txt
for fname in ~/GENR3/Imputed/1000G_PhaseIIIv5/chr*.info.gz
do
   zcat $fname | tail -n+2 -q >> sourcefiles/rsq.txt
   echo $'\n' >> sourcefiles/rsq.txt
done
sed -i '/^$/d' sourcefiles/rsq.txt
gzip sourcefiles/rsq.txt

# Merge chromosome-specific result files into a genome-wide file
# Main analysis ##################################################################
gzip -r results/mainanalysis/
zcat results/mainanalysis/chr22.SingleScore.assoc.gz | head -n1 > results/mainanalysis_autosomal.csv
for fname in results/mainanalysis/chr*.SingleScore.assoc.gz
do
   zcat $fname | tail -n+2 >> results/mainanalysis_autosomal.csv
done
gzip results/mainanalysis_autosomal.csv

# Secondary analysis: GIRLS ONLY ################################################ 
gzip -r results/girlsonly/
zcat results/girlsonly/chr22.SingleScore.assoc.gz | head -n1 > results/supp_girls_autosomal.csv
for fname in results/girlsonly/chr*.SingleScore.assoc.gz
do
   zcat $fname | tail -n+2 >> results/supp_girls_autosomal.csv
done
gzip results/supp_girls_autosomal.csv

# Secondary analysis: BOYS ONLY #################################################
gzip -r results/boysonly/
zcat results/boysonly/chr22.SingleScore.assoc.gz | head -n1 > results/supp_boys_autosomal.csv
for fname in results/boysonly/chr*.SingleScore.assoc.gz
do
   zcat $fname | tail -n+2 >> results/supp_boys_autosomal.csv
done
gzip results/supp_boys_autosomal.csv
