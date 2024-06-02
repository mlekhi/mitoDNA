# Analysis of genomic and epigenomic data from CLSA datasets

# Preprocesses the phenotype data
# Subsets the data to keep only rows with specified/significant haplotypes
# Convert certain variables to appropriate classes such as factors or numeric
# Integrates methylation data with phenotype and haplotype information
# Conducts statistical analysis using linear mixed effects regression models
# running lm for every single probe - lm(dependent variable, y ~ independent variable, covariates)
# Linear model fitting & p value calculation
# Outputs results (dataframe containing CpGs, coefficient estimates, standard errors, p-values) for further interpretation

# Dependencies: R programming language, data.table library

######################################
######   Setup + Data Loading  #######
######################################

# Set the working directory
setwd("/home/mlekhi/projects/def-ccastel/SharedResources/CLSA_private/results")
suppressMessages(library(data.table))

# Read in methylation data
bVals <- as.data.frame(fread("epigenetics/maya_bVals_rows.csv"))
dim(bVals) # 748702  1479
class(bVals)

# Read in phenotype data
pdata <- read.csv("epigenetics/pdataSVs.csv")

#################################
######  Data Preprocessing  #####
#################################

# Collapse smoking groups 2 and 3 together
pdata$smoke[which(pdata$smoke %in% c(2,3))] <- 2

# Remove missing smoking values and NA ly values
pdata <- pdata[-which(pdata$smoke == 99),]  
pdata <- pdata[-which(is.na(pdata$ly)),]

############################################

# Read haplotype data in
haplotype_data <- read.table("genomics/haplotypes.txt", header=TRUE)
head(haplotype_data)
dim(haplotype_data)

# Merge haplotype data with CLSA baseline data
clsa_baseline <- read.csv("baseline/b_epi_responses.csv")
# Merged based on GWAS linking key
matched_data <- merge(haplotype_data, clsa_baseline, by.x = "ADM_GWAS_COM", by.y = "ADM_GWAS3_COM")
dim(matched_data) # 1445   73

# Combine haplotype_data with phenotype data
pdata_haplo <- merge(matched_data, pdata, by.x = "ADM_EPIGEN2_COM")
dim(pdata_haplo) # 1301  110

# Define the columns to keep when trimming the dataframe for the analysis
desired_columns <- c("ADM_EPIGEN2_COM", "ADM_GWAS_COM", "haplotype", "age", "sex", "haplogroup", "smoke",
                     "ethnicity", "wbc", "ly", "mo", "gr", "row", "col", "SV1", "SV2", "SV3", "SV4", 
                     "SV5", "SV6", "SV7", "SV8", "SV9", "SV10", "entity_id", "IDAT_FILE_NAME_Grn_Cy3", 
                     "IDAT_FILE_NAME_Red_Cy5", "BS_Conv_Batch")

# Trim the pdata_haplo dataframe to keep only desired columns
trimmed_data <- pdata_haplo[, desired_columns]
dim(trimmed_data) # 1301  110

# Remove rows with missing values
pdata_haplo <- na.omit(trimmed_data)
dim(pdata_haplo) # 1301  28

#############################################

# Process bVals dataframe

# Remove "Red" suffix in pdata_haplo to create basename
pdata_haplo$Basename <- gsub("_Red.idat", "", pdata$IDAT_FILE_NAME_Red_Cy5)

# Remove "X" prefix in column names of bVals to create basename
colnames(bVals) <- sub("^X", "", colnames(bVals))

# Set column names as CpG sites
bVals <- t(bVals)
colnames(bVals) <- bVals[1, ]
bVals <- bVals[-1, ]

# Merging pdata haplo and bVals based on idat identifier 
pdata_haplo <- merge(pdata_haplo, bVals, by.x="Basename", by.y="row.names")
dim(pdata_haplo)

#######################################
######  Subsetting Selection ##########
#######################################

# Define the list of haplotypes to keep
haplotypes_to_keep <- c("110000", "100000", "000000", "000010", "000100", "000101", "001000", "010000")

# Subset the pdata_haplo dataframe to keep only rows with specified haplotypes
pdata_haplo <- pdata_haplo[pdata_haplo$haplotype %in% haplotypes_to_keep, ]
table(pdata_haplo$haplotype)

###############################################

# Fix classes of variables
haplotype <- as.factor(pdata_haplo$haplotype)
age <- as.numeric(pdata_haplo$age)
sex <- as.factor(pdata_haplo$sex)
smoke <- as.factor(pdata_haplo$smoke) 
BS_Conv_Batch <- as.factor(pdata_haplo$BS_Conv_Batch) 
row <- as.factor(pdata_haplo$row)
ly <- as.numeric(pdata_haplo$ly)
mo <- as.numeric(pdata_haplo$mo)
gr <- as.numeric(pdata_haplo$gr)
SV1 <- as.numeric(pdata_haplo$SV1)
SV2 <- as.numeric(pdata_haplo$SV2)
SV3 <- as.numeric(pdata_haplo$SV3)
SV4 <- as.numeric(pdata_haplo$SV4)
SV5 <- as.numeric(pdata_haplo$SV5)
SV6 <- as.numeric(pdata_haplo$SV6)
SV7 <- as.numeric(pdata_haplo$SV7)
SV8 <- as.numeric(pdata_haplo$SV8)
SV9 <- as.numeric(pdata_haplo$SV9)
SV10 <- as.numeric(pdata_haplo$SV10)

################################
###   Matching Participants  ###
################################

# Keep matched participants in bVals
common_row_names <- intersect(rownames(bVals), pdata_haplo$Basename)
bVals <- bVals[common_row_names, ]
bVals <- t(bVals)
dim(bVals) # 748702   1301

# Defining covariates
covariates <- c("age", "sex", "smoke", "ly", "mo", "gr", "SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7", "SV8", "SV9", "SV10", "row", "BS_Conv_Batch")

##################################
######  Statistical Analysis  ####
##################################

print("--- running lm")

## Define the lm function
# Define a function named runlm which takes a data frame as input
runlm <- function(thisdat) {
   # Evaluate the expression within the function
   lm1.out <- eval(parse(text=expression))
   # Get the summary of the linear model
   smodel = summary(lm1.out)
   return(smodel)
}

# Define a function named lmp which calculates the p-value for a given model object
lmp <- function (modelobject) {
    # Extract the F-statistic from the model object
    f <- modelobject$fstatistic
    # Calculate the p-value using the cumulative distribution function of the F-distribution
    p <- pf(f[1], f[2], f[3], lower.tail=F)
    # Remove attributes from the p-value object
    attributes(p) <- NULL
    return(p)
}

# Define variable names for matrices
varnames <- c("haplotype")

# Create matrices to store results
Bmat <- Rsmat <- Pmat <- matrix(NA, nrow=nrow(bVals), ncol=1)
rownames(Bmat) <- rownames(Rsmat) <- rownames(Pmat) <- rownames(bVals)

colnames(Bmat) <- paste("Estimate", varnames, sep=".")
colnames(Rsmat) <- paste("R-squared", varnames, sep=".")
colnames(Pmat) <- paste("p-value", varnames, sep=".")

# Loop through each row of the bVals dataframe
for (i in 1:nrow(bVals)) {
  # Print progress message for every 100,000th iteration
  if (i %% 100000 == 0) {cat(paste("On probe ", i, "\n", sep=""))}
  # Extract the expression data for the current row
  thisExpr <- as.numeric(bVals[i,])
  class(thisExpr)
  class(haplotype)
  # Construct the linear model expression with covariates
  expression <- "lm(thisExpr ~ haplotype 
    + age + sex + smoke + ly + mo + gr + 
    SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + 
    SV8 + SV9 + SV10 + row + BS_Conv_Batch)"

  # Create a design matrix including expression and covariates
  designmatrix <- data.frame(thisExpr, haplotype, age, sex, smoke, ly, mo, gr, SV1, SV2, SV3, SV4, SV5, SV6, SV7, SV8, SV9, SV10, row, BS_Conv_Batch)
  # Running linear model
  lm1.out <- try(runlm(designmatrix), silent=F)

  # Check if the linear model was successfully computed
  if (substr(lm1.out[1], 1, 5) != "Error") {
    # Extract coefficients from the linear model summary
    tabOut <- lm1.out$coefficients
    # Store the coefficient estimates, standard errors, and p-values
    Bmat[i,] <- tabOut[2,"Estimate"]
    Rsmat[i,] <- tabOut[2,"Std. Error"]
    Pmat[i,] <- tabOut[2,"Pr(>|t|)"]
  } else {
      # Printng NA erorr message
      cat('Setting P-value=NA, Beta value=NA, and SE=NA\n')
      # Set values to NA in case of error
      Bmat[i,] <- Rsmat[i,] <- Pmat[i,] <- NA
  }
}

# Combine coefficient estimates, standard errors, and p-values into a dataframe
lm_res <- cbind(Bmat, Rsmat, Pmat)
lm_res <- as.data.frame(lm_res)

###########################
#####   Saving Data   #####
###########################

# Writing results to a CSV file
write.csv(lm_res, "epigenetics/lm_res_haplotype_noresids_v5.csv", quote=F)
print("--- DONE")

# Save R environment
save.image("/home/mlekhi/projects/def-ccastel/mlekhi/noresids_v5.RData")
