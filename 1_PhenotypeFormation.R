# Generation R analyses
# Within GENA, navigate to your project folder and open this script in Rstudio 
# (type rstudio in the terminal to launch) and run interactively 
# SCRIPT 1, what this script does:
#   - Read in phenotype, sex, age, and genetic principal component covariate data
#   - Create the correct variables, scale outcomes, and take out non-complete cases, non-Europeans and siblings
#   - Merge data into phenotype and covariates files
#   - Check sample descriptives

# Load packages and set directory
library(foreign) # To read spss files
library(psych) # For descriptives

# Read in phenotype data created in SCRPT 0. This contains the following variables
#   - IDC, 
#   - MOTHER: mother id (for sibling exclusion) you can find in the general data files
#   - GWAv3European or GWAv4European: indicating which children have been genotyped data, 
#     you can find the variable in the PCA_Selection_GWAv3*.sav or PCA_Selection_GWAv4*.sav files, 
#     note that here we are including European-ancestry participants only.
#   - phenotype: the score or variable of interest
#   - GENDER: child sex at birth (1 for boys and 2 for girls) you can also find in the general data files
#   - all other covariates: e.g. we will use child age and the first 10 genetic ancestry components 
#     (only European-ancestry participants) you can find in the PCA_Selection_GWAv3*.sav or PCA_Selection_GWAv4*.sav files
pheno.data <- read.csv("./sourcefiles/phenotype_data.csv")

# Get some descriptive to make sure nothing funky is going on
print(psych::describe(pheno.data[, c('phenotype', 'age')]))

# Create variables needed for RVtest (format mostly like in Plink: http://zzz.bwh.harvard.edu/plink/data.shtml)
pheno.data$fid <- pheno.data$IDC # Create family ID (siblings will be deleted -> family ID is identical to IDC)
pheno.data$iid <- pheno.data$IDC # Participant ID
pheno.data$fatid <- 0 # No parental data, so put 0
pheno.data$matid <- 0 # No maternal data, so put 0
pheno.data$sex <- as.numeric(as.factor(pheno.data$GENDER)) - 1 # Gender variable with 0 for boys and 1 for girls

#-------------------------------------------------------------------------------
# Exclusion steps: genotype data available: 
genot <- subset(pheno.data, !is.na(GWAv3European)); message("Original sample = ", nrow(genot)) 

# Check how many kids are excluded because of missing phenotype or covariate data
# Note: all kids' sex is known, so no exclusions needed
out <- genot[!is.na(genot[, 'phenotype']), ]
message(nrow(out), ' (', nrow(genot),' - ',nrow(out),' = ', nrow(genot)-nrow(out),' kids without outcome data')
age <- out[!is.na(out[, 'age']), ]
message(nrow(age), ' (', nrow(out),' - ',nrow(age),' = ', nrow(out)-nrow(age),' kids with unknown age')

#-------------------------------------------------------------------------------
# Only keep relevant variables and exclude non-complete cases.  
keep <- c("fid","iid","fatid","matid","sex",'MOTHER', paste0("C", 1:10), 'age' , 'phenotype')
pheno.data <- pheno.data[complete.cases(pheno.data[c(keep)]), c(keep)]

# remove related individuals
rm_siblings <- function(dataset) {
  sib_ids <- dataset[duplicated(dataset$MOTHER) | duplicated(dataset$MOTHER, fromLast = T), 'iid']
  outdata <- dataset[!dataset$iid %in% sib_ids, ]
  message(nrow(outdata), ' (', nrow(dataset),' - ',length(sib_ids),' siblings or twins)')
  return(outdata)
}
pheno.data  <- rm_siblings(pheno.data)

# ------------------------------------------------------------------------------
# Sex count 
message(nrow(subset(pheno.data, sex == 0)), " boys + ", nrow(subset(pheno.data, sex == 1)), " girls")

#-------------------------------------------------------------------------------
# OPTIONAL: Scale the phenotype (note this should be done only now, that we have the analytical sample, 
# so mean is precisely 0 and SD 1)
pheno.data[,'phenotype_z'] <- scale(pheno.data[, 'phenotype'])
pheno.data[,'age_z'] <- scale(pheno.data[, 'age'])
print(psych::describe(pheno.data[, c(ncol(pheno.data)-1,ncol(pheno.data))]))

# ------------------------------------------------------------------------------

write_output <- function(dataset){
  # Write outcome and covariate files
  write.csv(dataset, file='./sourcefiles/phenotype.csv')
  write.table(dataset[c(keep[1:5], 'phenotype', 'phenotype_z')], 
              file = './sourcefiles/phenotype.ped', quote = F, row.names = F)
  write.table(subset(dataset, sex == 0, select = c(keep[1:5], 'phenotype', 'phenotype_z')), 
              file = './sourcefiles/phenotype_boys.ped', quote = F, row.names = F)
  write.table(subset(dataset, sex == 1, select = c(keep[1:5], 'phenotype', 'phenotype_z')), 
              file = './sourcefiles/phenotype_girls.ped', quote = F, row.names = F)
  write.table(dataset[c(keep[1:5], 'age', paste0("C", 1:10))], 
              file = './sourcefiles/covariates.ped', quote = F, row.names = F)
}
write_output(pheno.data)

##########################
# For details of QC, see:  
# 1) Medina-Gomez C, Felix JF, Estrada K, Peters MJ, Herrera L, Kruithof CJ et al. Challenges in conducting 
# genome-wide association studies in highly admixed multi-ethnic populations: the Generation R Study. Eur J Epidemiol 2015; 30: 317330.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4385148/ 
# 2) Lamballais S et al. Genetic scores for adult subcortical volumes associate with subcortical volumes during infancy and childhood.
# Human Brain Mapping 2021; 42: 1583-93. https://doi.org/10.1002/hbm.25292

descr <- function(dataset) {
  # get most descriptives
  t <- as.data.frame(psych::describe(dataset[, 'phenotype']))
  t['girls',] <- as.data.frame(psych::describe(subset(dataset, sex==1)['phenotype']))
  t['boys',]  <- as.data.frame(psych::describe(subset(dataset, sex==0)['phenotype']))
  t['age',]   <- as.data.frame(psych::describe(dataset[, 'age']))
  # get IQR
  iqr <- data.frame("IQR" = IQR(dataset[, 'phenotype']))
  iqr['girls',] <- IQR(subset(dataset, sex==1)[, 'phenotype'])
  iqr['boys',]  <- IQR(subset(dataset, sex==0)[, 'phenotype'])
  iqr['age',]   <- IQR(dataset[, 'age'])
  # bind all together
  t1 <- cbind(t[c('n','mean','sd','median')], iqr, t[c('min','max','skew','kurtosis','se')])
  t1 = round(t1, 5) # show 5 decimal places
  View(t1)
}
descr(pheno.data)

# Histograms -------------------------------------------------------------------

title = "==> Insert title"

pdf("./results/OutcomeHistograms.pdf")

hist(pheno.data$phenotype, main = title, xlab = "Unscaled phenotype")
hist(pheno.data[pheno.data$sex == 1, 'phenotype'], 
     main = paste(title, '(girls)'), xlab = "Unscaled phenotype")
hist(pheno.data[pheno.data$sex == 0, 'phenotype'], 
     main = paste(title, '(boys)'), xlab = "Unscaled phenotype")

hist(pheno.data$fm_score_6m_z, main = title, xlab = "Scaled phenotype")
hist(pheno.data[pheno.data$sex == 1, 'phenotype_z'], 
     main = paste(title, '(girls)'), xlab = "Scaled phenotype")
hist(pheno.data[pheno.data$sex == 0, 'phenotype_z'], 
     main = paste(title, '(boys)'), xlab = "Scaled phenotype")

dev.off()

