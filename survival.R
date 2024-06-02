# Analyzing significant CpGs from reggression using gene ontology (GO) enrichment analysis

# Dependencies: R programming language, data.table library, BiocManager, missMethyl BioConductor package

######################################
######   Setup + Data Loading  #######
######################################

# Set the working directory
setwd("/home/mlekhi/projects/def-ccastel/SharedResources/CLSA_private/results")

# Check if the BiocManager package is installed and install it if not
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install the missMethyl package using BiocManager
BiocManager::install("missMethyl")

# Read the regression results into a data frame
suppressMessages(library(data.table))
lm <- as.data.frame(fread("epigenetics/lm_res_haplotype_noresids_v5.csv"))

# Rename the fourth column for ease of access
colnames(lm)[4] <- "pHaplotype"

# Extract all CpGs from the linear regression results
allcgps = lm$V1
length(allcgps)
class(allcgps)

# Extract only significant CpGs based on the significance level
sigcpgs = subset(lm, lm$pHaplotype < 1E-07)
dim(sigcpgs)
# Gathering IDs of signficant CpGs
sigcpgs <- sigcpgs$V1
length(sigcpgs)
class(sigcpgs)


####################################
######  Running GO Analysis    #####
####################################

# Run GO analysis using the gometh function
# Setting significant CpGs to list from prev step
# Submitting all CpGs from the regression results
gst.go.all <- gometh(sig.cpg = sigcpgs, all.cpg = allcgps, 
                     collection = "GO", prior.prob = TRUE, array.type = c("EPIC"),
                     genomic.features = c("ALL"))

# Order the results by FDR to identify significant results & for better interpretation
gst.go.all.orderedfdr <- gst.go.all[order(gst.go.all$FDR), ]