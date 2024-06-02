# Filtering and processing epigenetic data obtained from microarray experiments 
# 1. Quality Control - removing poor quality samples (duplicates - maintaining those with higher call rate, filter for high call rate - removing high detection p value)
# 2. Normalization - background correction
#    noob normalization was selected for its ability to remove variation between samples in the dataset (quantile & noob were tested)
# 3. Probe-level Filtering - removing cross reactive probes, removing sex probes, removing SNPs at CpG site
# 4. Methylation Beta Values Calculation - storing computed beta values in in data frame with CpG site info

# Dependencies: R programming language, minfi package, Illumina Human Methylation EPIC array annotation package, data.table library

setwd("/home/mlekhi/projects/def-ccastel/SharedResources/CLSA_private/data/epigenetics")

# Load packages
suppressMessages(library(minfi))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

############################
### Load data into minfi ###
############################

# Read in the sample sheet for the experiment
dataDirectory <-  "/home/mlekhi/projects/def-ccastel/SharedResources/CLSA_private/data/epigenetics"
pdata <- read.metharray.sheet(dataDirectory, pattern="CLSA_metadata_FINAL_Apr2021-deidentified.csv")
dim(pdata) # 1478   17

# Basename is used to feed minfi the prefix of the idat
# Getting rid of the words red.idat
pdata$Basename <- gsub("_Red.idat", "", pdata$IDAT_FILE_NAME_Red_Cy5)

# Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(base = dataDirectory, targets=pdata)
sampleNames(rgSet) <- pdata$Basename

# Calculate the detection p-values
detP <- suppressMessages(detectionP(rgSet))


############################
##    Quality Control   ####
############################

# Track how many samples were removed due to poor quality
samples_before <- dim(rgSet)[2]

# Remove poor quality samples with high detection p value
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

# Print out number of samples removed due to poor quality
samples_removed <- samples_before - dim(detP)[2]
message("----- ", samples_removed, " sample(s) removed due to poor quality") #----- 0 sample(s) removed due to poor quality


############################
####    Normalization   ####
############################

# Normalize the data using minfi's preprocess noob (normal-exponential out-of-band background correction)
mSetSq <- suppressMessages(preprocessNoob(rgSet))
mSetSq <- mapToGenome(mSetSq)
mSetSq <- ratioConvert(mSetSq)

#################################
#####  Probe-level Filtering ####
#################################

# Ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# Remove any probes that have failed in one or more samples
probes_before <- dim(mSetSq)[1]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

# Print out number of probes removed for failing in one or more samples
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for failing in one or more samples") #----- 35533 probe(s) removed for failing in one or more samples
probes_before <- dim(mSetSqFlt)[1]

############################
# Remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

# Print out number of probes removed for having SNPs at CpG site
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for having SNPs at CpG site") #----- 25504 probe(s) removed for having SNPs at CpG site
probes_before <- dim(mSetSqFlt)[1]

############################
# Exclude cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "PidsleyCrossReactiveEPIC.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]

# Print out number of probes removed for being cross reactive
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for being cross reactive") #----- 24846 probe(s) removed for being cross reactive
probes_before <- dim(mSetSqFlt)[1]


############################
# Remove Sex Probes
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")),]
keep <- !(featureNames(mSetSqFlt) %in% sexProbes$Name)
mSetSqFlt <- mSetSqFlt[keep,]

# Print out number of probes removed for being on sex chromosomes
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for being on sex chromosomes") #----- 16805 probe(s) removed for being on sex chromosomes

# Print out the number of probes remaining
message("----- ", dim(mSetSqFlt)[1], " probe(s) remaining for analysis") #----- 763171 probe(s) remaining for analysis

########################################################
########  Methylation Beta Values Calculation ##########
########################################################

# Calculate methylation beta values
bVals <- getBeta(mSetSqFlt)
head(bVals)
dim(bVals)

# writing results to a CSV
write.csv(bVals, file = "../../maya_bVals_rows.csv", quote = FALSE)

save.image("/home/mlekhi/projects/def-ccastel/mlekhi/QC.RData")