# Correlation Analysis with Haplotypes and Variables
# Constructed to examine the interplay between age, sex, smoking, white blood cell counts, row, batch, and the initial 10 SVs
# Purpose: explore relationships and potential associations between haplotypes and phenotypic traits/other covariates

# Dependencies: R programming language, data.table library


##########################
#### Data Preparation ####
##########################

# Load base dataset
load("/home/mlekhi/projects/def-ccastel/mlekhi/correlation_base_data.RData")

# Select desired columns for analysis
desired_columns <- c("haplotype", "age", "sex", "smoke",
                     "ly", "mo", "gr", "SV1", "SV2", "SV3", "SV4", 
                     "SV5", "SV6", "SV7", "SV8", "SV9", "SV10", "BS_Conv_Batch")

# Subsetting the data frame to retain only selected columns
pdata_haplo <- pdata_haplo[, desired_columns]

# Convert all columns to numeric
pdata_haplo$haplotype <- as.numeric(pdata_haplo$haplotype)

# Setting sex to 1 (male) and 2 (female)
pdata_haplo$sex <- ifelse(pdata_haplo$sex == "M", 1, ifelse(pdata_haplo$sex == "F", 2, NA))
pdata_haplo$sex <- as.numeric(pdata_haplo$sex)
pdata_haplo$age <- as.numeric(pdata_haplo$age)

pdata_haplo$smoke <- as.numeric(pdata_haplo$smoke) 
pdata_haplo$BS_Conv_Batch <- as.numeric(pdata_haplo$BS_Conv_Batch) 
pdata_haplo$ly <- as.numeric(pdata_haplo$ly)
pdata_haplo$mo <- as.numeric(pdata_haplo$mo)
pdata_haplo$gr <- as.numeric(pdata_haplo$gr)
pdata_haplo$SV1 <- as.numeric(pdata_haplo$SV1)
pdata_haplo$SV2 <- as.numeric(pdata_haplo$SV2)
pdata_haplo$SV3 <- as.numeric(pdata_haplo$SV3)
pdata_haplo$SV4 <- as.numeric(pdata_haplo$SV4)
pdata_haplo$SV5 <- as.numeric(pdata_haplo$SV5)
pdata_haplo$SV6 <- as.numeric(pdata_haplo$SV6)
pdata_haplo$SV7 <- as.numeric(pdata_haplo$SV7)
pdata_haplo$SV8 <- as.numeric(pdata_haplo$SV8)
pdata_haplo$SV9 <- as.numeric(pdata_haplo$SV9)
pdata_haplo$SV10 <- as.numeric(pdata_haplo$SV10)

###################################
###   Correlation Analysis    #####
###################################

# Setting new working directory
setwd("/home/mlekhi/projects/def-ccastel/mlekhi/correlation-matrices")

# Open a PDF file for storing correlation matrices
pdf("pcCorMatrix_noresids_v2.pdf")

# Define function to plot correlation matrix
twolines = function(x,y) {
  # using the haplotype column
  points(x,y,pch=20, col=pdata_haplo$haplotype)
  abline(lm(y~x),col="red")
}

#   x: Numeric vector for the x-axis
#   y: Numeric vector for the y-axis
#   digits: Number of digits to display for correlation coefficient (default is 2)
#   prefix: Prefix to be added before the correlation coefficient (default is "")
#   cex.cor: Size of the correlation coefficient text (default is calculated based on the length of the text)
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
#   x: Numeric vector for the x-axis
#   labels: Labels for the plot
mydiag.panel <- function( x,  labels, ...) {
  ll <- par("usr")
  rect(ll[1], ll[3], ll[2], ll[4], col="darkolivegreen1")
}

diag.labels=c("haplotype", "age", "sex", "smoke", "ly", "mo", "gr", "SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7", "SV8", "SV9", "SV10", "BS_Conv_Batch")

# Define plot formula
plot.formula=as.formula(~haplotype + age + sex + smoke + ly + mo + gr + SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10)

pairs(plot.formula, data=pdata_haplo, upper.panel=twolines, labels=diag.labels, diag.panel=mydiag.panel, lower.panel=panel.cor, label.pos=0.5, main="Correlation Between Variables")
dev.off()