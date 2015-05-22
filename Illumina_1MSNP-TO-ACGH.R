# Author Wibo Pipping 22-05-2015
# contact: wibo@thehyve.nl
# This R script describes how to transform the 1M Illumina SNP data for upload in aCGH format.
# Of the 9 input files only 4 will be used. calls.txt, snpsegments.txt, segments.txt and probes.txt
# Using the probes.txt file platforms are generated that can be uploaded to tranSMART.
# The calls.txt file contains call values for big loss, loss, gain or big gain.
# snpsegments and segments are the log ratio and B Allele frequency.
# NOTE: calls data can be loaded and used in some of the ACGH analysis, segment and snpsegment data cannot.

########################################################################################################################
# Count function needed to produce the platform regions

# getSubset selects a range of probes that calls 
getSubset <- function(df2=calls, df1=probes) {
  tempSubset <- df1[df1[,"Chromosome"]==df2[,"Chromosome"],]
  numOfProbes <- sum(df2[,"Start"] <= tempSubset$End & tempSubset$End <= df2[,"End"])
  return(numOfProbes)
}
########################################################################################################################
# Loading the data in R, variables needed
# probes
# segments
# snpsegments
# calls

probes <- read.csv("", header = T , sep ="\t", colClasses=c("character","numeric", "numeric", "numeric", "numeric"))
calls <- read.csv("", header = T , sep ="\t", colClasses=c("character","numeric", "numeric", "numeric"))
segments <- read.csv("", header = T , sep ="\t", colClasses=c("character","numeric", "numeric", "numeric"))
snpsegments <- read.csv("", header = T , sep ="\t", colClasses=c("character","numeric", "numeric", "numeric"))

# Variables needed to build the platform file
# CYTOBAND, GENE_SYMBOL, GENE_ID and ORGANISM
# empty values and Homo Sapiens as I do not have this information
cytoband <- ""
geneSym <- ""
geneId <- ""
organism <- 'Homo Sapiens'
########################################################################################################################
# build calls platform

# number of probes in a region. First specifies an empty vector with NA values for each row in the dataset
# next counts the amount of probes found in a region using the getSubSet() function defined above.
numProbes <- rep(NA, nrow(calls))
for(i in seq(1:nrow(calls))) {
  numProbes[i] <- getSubset(calls[i,])
}

# Set GPL_ID, which is the platform name stored in the database.
gplID <- rep("1M_SNP_Calls", nrow(calls))

# Generate region names based on the data 
regions <- rep(NA,nrow(calls))
for(x in seq(1:nrow(calls))) {
  regions[x] <- paste(calls[x,1],paste(calls[x,2],calls[x,3],sep="-"),sep=":")
}

#  Build platform file with 
calls_platform <- cbind(gplID,regions,substring(calls[,1],first=4),calls[,2:3],numProbes,cytoband,geneSym,geneId,organism)
colnames(calls_platform) <- c("GPL_ID","REGION_NAME", "CHROMOSOME","START_BP","END_BP","NUM_PROBES","CYTOBAND","GENE_SYMBOL","GENE_ID","ORGANISM")
head(calls_platform)

# Export data
write.table(calls_platform, file="calls_platform.txt", quote=F, sep="\t", row.names=F, append=F)
#write.table(calls_platform, file="~/transmart-data/samples/studies/Cell-line/Illumina_1M_SNP_Array_from_Nexus/calls_platform.txt", quote=F, sep="\t", row.names=F, append=F, na="")

########################################################################################################################
# build calls data file

# upload requires fields to have data, so empty values are replaced with 0
array.chip <- 0
array.region <- 0
array.loss <- 0
array.norm <- 0
array.gain <- 0
array.amp <- 0

# build the calls data file
aCGH_calls <- cbind(regions, array.chip, array.region, calls[,4],array.loss,array.norm,array.gain,array.amp)
colnames(aCGH_calls) <- c("region_id","array.chip", "array.segm", "array.flag", "array.loss", "array.norm", "array.gain", "array.amp")

# Write calls data file
write.table(aCGH_calls, file="aCGH_calls.txt", quote=F, sep="\t", row.names=F, append=F, na="")
#write.table(aCGH_calls, file="~/transmart-data/samples/studies/Cell-line/Illumina_1M_SNP_Array_from_Nexus/pc346_calls.txt", quote=F, sep="\t", row.names=F, append=F, na="")

########################################################################################################################
# B Allele frequency
# build platform file, Left to right:
# "GPL_ID","REGION_NAME", "CHROMOSOME","START_BP","END_BP","NUM_PROBES","CYTOBAND","GENE_SYMBOL","GENE_ID","ORGANISM"

#GPL_ID
gplID <- rep("1M_SNP", nrow(segments))

# REGION_NAME
regions <- rep(NA, nrow(segments))
for(x in seq(1:nrow(segments))) {
regions[x] <- paste(segments[x,1],paste(segments[x,2],segments[x,3],sep="-"),sep=":")
}

# CHROMOSOME, START_BP, END_BP is taken from the input dataframe

# NUM_PROBES, uses the getSubset() function 
numProbes <- rep(NA, nrow(segments))




# Assemble Platform file ready for export
platform <- cbind(gplID,regions,substring(segments[,1],first=4),segments[,2:3],numProbes,cytoband,geneSym,geneId,organism)
colnames(platform) <- c("GPL_ID","REGION_NAME", "CHROMOSOME","START_BP","END_BP","NUM_PROBES","CYTOBAND","GENE_SYMBOL","GENE_ID","ORGANISM")
head(platform)

# export
write.table(platform, file="~/transmart-data/samples/studies/Cell-line/Illumina_1M_SNP_Array_from_Nexus/1M_SNP_PC346c_platform.txt", quote=F, sep="\t", row.names=F, append=F)

########################################################################################################################
# Rebuild Segments and SNPsegments to fit in aCGH data format.

# Empty placeholders
chip.flag <- 0
chip.loss <- 0
chip.norm <- 0
chip.gain <- 0
chip.amp <- 0

aCGH_Data <- cbind(regions, as.numeric(segments[,4]), as.numeric(snpsegments[,4]),chip.flag, chip.loss, chip.norm, chip.gain, chip.amp)
colnames(aCGH_Data) <- c("region_id","array.chip", "array.segm", "array.flag", "array.loss", "array.norm", "array.gain", "array.amp")
head(aCGH_Data)

# export
write.table(aCGH_Data, file="~/transmart-data/samples/studies/Cell-line/Illumina_1M_SNP_Array_from_Nexus/pc346_acgh_data.txt", quote=F, sep="\t", row.names=F, append=F, na="")
