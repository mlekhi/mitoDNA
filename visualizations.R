# Scripts for the analysis of epigenetic data 
# Visualization types: QQ plot, manhattan plot, volcano plot

# Note: adjust input and output file paths as needed based on your data
# Dependencies: R programming language, data.table library, QCEWAS package, ggplot2 package, IlluminaHumanMethylationEPICanno.ilm10b4.hg19

#############################
####  Loading in Data   #####
#############################

# Set working directory to the results folder
setwd("/home/mlekhi/projects/def-ccastel/SharedResources/CLSA_private/results")

# Load required libraries with messages suppressed
suppressMessages(library(data.table))
suppressMessages(library(QCEWAS))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(ggplot2))

lm <- as.data.frame(fread("epigenetics/lm_res_haplotype_noresids_v5.csv"))
# fixing col name p.data.haplotype
# make sure to use raw, not adjusted for volcano plot

# Rename haplotype column for ease of access
colnames(lm)[4] <- "pHaplotype"

# Read EPIC array annotation data
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- as.data.frame(ann)

# Merge linear regression results with annotation data
lm_res <- merge(lm, ann, by.x="V1", by.y="row.names")

# Multiple test correction using FDR
lm_res$adj.P.Val <- p.adjust(lm_res$pHaplotype, method = "BH", n = length(lm_res$pHaplotype))

#calculate lambda value
# note: you may need to rename the column
pvals <- lm_res$pHaplotype
# pvals <- lm$p.value.haplotype
P_lambda(pvals) # 1.547072; noresids 1.356139; noresids_v2: 0.8133869; noresids_v3: 0.7668753; noresids_v4: 1.1404, noresids_v5: 1.148973


##############
### QQplot ###
##############

observed <- sort(pvals)
observed2 <- c(length(pvals))
observed2null <- -(log10(observed2 / (length(observed2)+1)))
pvals <- c(pvals, observed2null)
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
m="qqplot"

# Creating a PDF file to save the QQ plot
pdf(file="epigenetics/qqplot_maya_noresids_v5.pdf", width=10)

# Plotting the QQ plot with red diagonal line
plot(c(0,20), c(0,20), col="red", lwd=4, type="l", xlab="Expected (-logP)",
     ylab="Observed (-logP)", xlim=c(0,16), ylim=c(0,16), las=1, xaxs="i", yaxs="i", bty="l", main=m)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")

dev.off()


######################
### Manhattan plot ###
######################

# Rename colnames of lm_res
colnames(lm_res)[1] <- "V1"
colnames(lm_res)[2] <- "estimate"
colnames(lm_res)[4] <- "P.Value"
colnames(lm_res)[6] <- "start"

# Create a new column - Rounding the "start" column to the nearest integer
lm_res$position<-round(lm_res$start,digits=0)

# Convert chr format to format without the leading 'chr'
lm_res$chr <- gsub("^.{0,3}", "", lm_res$chr)

# Define a function named "manhattan" to create a Manhattan plot
manhattan<-function(data, sig=NULL){
  # Subset data to remove chromosomes "0", "X", and "Y"
  data<-subset(data, (data$chr!="0"))
  # Remove rows where "start" column has "NA" values
  data<-subset(data, (data$start!="NA"))
  data<-subset(data, (data$chr!="X"))
  data<-subset(data, (data$chr!="Y"))

  # Remove "chr" prefix and convert "chr" column to numeric
  data$chr <- gsub("chr","",data$chr)
  data$chr <- as.numeric(data$chr)

  # Order the data by chromosome number
  data2 <- data[order(data$chr),]
  data=data2

  # Calculate the minimum y-axis value for the Manhattan plot
  if (is.null(sig)==TRUE){
    ymin1=round(max(-log10(data$P.Value))+1)
  }
  else{
    ymin1=34
  }

  # Starting a PDF for the manhattan plot
  title=c()
  pdf("epigenetics/manhattan_plot_noresids_v5.pdf",width = 45,height = 18)

  # Setting chromosome count
  chr <- c(1:22)

  #Summary statistics
  print(table(data$chr))
  par(mar=c(5,5,2,2))
  phy.max<-tapply(data$start, data$chr,max,na.rm=T)
  cumlen=0
  
  # Iterate over chromosomes to calculate cumulative position
  for(i in chr){
    cat(paste("Now working on chromosome ",i,"\r"))
    data[data$chr==i,"loc"]<-data[data$chr==i,"position"]+cumlen
    cumlen<-cumlen+phy.max[i]
  }
  phy.med<-tapply(data$loc,data$chr,median,na.rm=T)[chr]
  print(phy.med)
  data$mlgpval<- -log(data[,"P.Value"], base=10)
  
  # creating blank plot
  plot(data[,"loc"],data[,"mlgpval"],type="n",yaxt="n",xaxt="n",
       xlab="Chromosome",
       ylab=expression(-log[10]*P),main=title,
       xlim=c(0,max(data$loc,na.rm=T)),cex.lab=2.5,ylim=c(0,ymin1))
  col=ifelse(data$mlgpval< 10, (rep(c("#4B4059","#57486F"),13)), "#BD6B73")
  
  # Adding axis ticks
  axis(side=2, at=seq(from=0,to=ymin1,by=2), labels=seq(from=0,to=ymin1,by=2),
       tick=T,cex.axis=2,las=1)
  axis(side=1, at=phy.med[c(1:22)], labels=chr[c(1:22)],
       tick=T,cex.axis=2,las=1)

  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white")

  # Add points to the plot for each chromosome
  for(i in chr){
    cat(paste("Now working on chromosome ",i,"\r"))
    points(data[data$chr==i,"loc"],data[data$chr==i,"mlgpval"],
           col=ifelse(data[data$chr==i,"mlgpval"] > sig, "#BD6B73", col[i]),pch=20,cex=2)
  }
  
  # Add a dotted line for significance threshold if provided
  if (is.null(sig)==FALSE){
    abline(h=sig,lty="dotted",lwd=5,col="#57486F")
  }
  
  dev.off()
}

# Running manhattan function
manhattan(data = lm_res, sig=(-log10(1e-06)))


####################
### Volcano Plot ###
####################

# Renaming column
colnames(lm)[4] <- "pHaplotype"

png("epigenetics/volcano_plot_noresids_v5.png", width=250, height=250, units='mm', res = 600, pointsize=80)
ggplot() +
  geom_point(data = lm_res, aes(x = estimate, y = -log10(P.Value)), color = "#4B4059") +
  geom_point(data = lm_res[which(lm_res$pHaplotype < 1e-06),], aes(x = estimate, y = -log10(P.Value)), color = "#BD6B73") +
  geom_hline(yintercept = -log10(1e-06), linetype = "dotted", color = "#57486F", size = 0.75) +
  xlim(-0.08, 0.08) + ylim(0, 16) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 17), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=17, vjust = 0.5), 
        axis.text.y = element_text(size=17)) +
  labs(x = "Beta Estimate", y = expression(-log[10]*P)) 
dev.off()
    

##############################
### Alternate Volcano Plot ###
##############################

lm_res$col <- ifelse(lm_res$pHaplotype < 1E-06, "red", "black")
DMS2 <- subset(lm_res, lm_res$col == "red")

volcano <- function(DMS, filename, sig = NULL) {
  pdf(paste("volcano_", filename, ".pdf", sep = ""), width = 16, height = 16)
  
  plot(-log10(DMS$pHaplotype) ~ DMS$Estimate.haplotype,
       xlab = "Beta Estimate",
       ylab = expression(-log[10]*P),
       main = "Volcano Plot", ylim = c(0, 10), xlim = c(-0.5, 0.5), pch = 20, col = DMS$col)
  
  par(new = TRUE)
  
  plot(-log10(DMS2$pHaplotype) ~ DMS2$Estimate.haplotype, axes = FALSE,
       ylab = expression(-log[10]*P), xlab = "",
       ylim = c(0, 10), xlim = c(-0.5, 0.5), pch = 20, col = DMS2$col)
  
  if (!is.null(sig)) {
    abline(h = sig, lty = "dotted", lwd = 2, col = "chartreuse4")
  }
  
  dev.off()
}

# Assuming DMS is the data frame from lm_res
volcano(lm_res, filename = "noresids_v2", sig = 7)