################################################################
################################################################
# Script: ProgR_1.1_SNP.R
# Author: Ingrid Lonnstedt
# Date:  15/02/2013
# R version: R 2.15.1
# Details: SNP analysis
################################################################
################################################################



################################################################
################################################################
#
# File paths, functions and date
#
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.raw = paste(getwd(), "/rawData/responsify/HG-U133_Plus_2/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/HG-U133_Plus_2/", sep='')
prefix.out = paste(getwd(), "/output/Expression/", sep='')

#prefix.prog = paste(getwd(), "/programs/", sep='')
#source(paste(prefix.prog, 'ProgR_0.1_Functions.R', sep=''))
date = format(Sys.Date())
library(aroma.affymetrix)

### Fileinfo: only the files with pass == TRUE to be used
###
filenames = list.files(substr(prefix.raw, 1, nchar(prefix.raw, type = 'c')-1))
pass = ifelse(substr(filenames, 1, 1) == 'B', T, F)


################################################################
################################################################
#
# Aroma Read and Preprocess files
#
################################################################
################################################################



verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

#Read cdf file
cdf <- AffymetrixCdfFile$byChipType("HG-U133_Plus_2");
print(cdf);

AffymetrixCdfFile:
  Path: annotationData/chipTypes/HG-U133_Plus_2
Filename: HG-U133_Plus_2.cdf
File size: 25.04 MB (26251783 bytes)
Chip type: HG-U133_Plus_2
RAM: 0.00MB
File format: v4 (binary; XDA)
Dimension: 1164x1164
Number of cells: 1354896
Number of units: 54675
Cells per unit: 24.78
Number of QC units: 9



#Read CEL files

# Raw CEL files using the new combined cdf (i.e. includes controls)
csR <- AffymetrixCelSet$byName('responsify', cdf=cdf); 
print(csR);

AffymetrixCelSet:
  Name: responsify
Tags: 
  Path: rawData/responsify/HG-U133_Plus_2
Platform: Affymetrix
Chip type: HG-U133_Plus_2
Number of arrays: 97
Names: 01-9528_(HG-U133_Plus_2), 03-12391_(HG-U133_Plus_2), 03-9183_(HG-U133_Plus_2), ..., BO5649_(HG-U133_Plus_2) [97]
Time period: 2011-10-04 11:53:19 -- 2012-10-10 15:26:51
Total file size: 1255.13MB
RAM: 0.12MB

#Proceed only with arrays that have info$pass == TRUE
csR = extract(csR, (1:length(filenames))[pass])
print(csR)

AffymetrixCelSet:
  Name: responsify
Tags: 
  Path: rawData/responsify/HG-U133_Plus_2
Platform: Affymetrix
Chip type: HG-U133_Plus_2
Number of arrays: 86
Names: BO2098_(HG-U133_Plus_2), BO2196_(HG-U133_Plus_2), BO2229_(HG-U133_Plus_2), ..., BO5649_(HG-U133_Plus_2) [86]
Time period: 2011-10-04 11:53:19 -- 2012-10-10 15:26:51
Total file size: 1112.81MB
RAM: 0.10MB

ae <- ArrayExplorer(csR)
setColorMaps(ae, "sqrt,yellow")
process(ae, verbose=TRUE)
print(ae)


# background corrections RMA
bc <- RmaBackgroundCorrection(csR, tag="BC");
csBC <- process(bc,verbose=verbose);
print(csBC);

AffymetrixCelSet:
  Name: responsify
Tags: BC
Path: probeData/responsify,BC/HG-U133_Plus_2
Platform: Affymetrix
Chip type: HG-U133_Plus_2
Number of arrays: 86
Names: BO2098_(HG-U133_Plus_2), BO2196_(HG-U133_Plus_2), BO2229_(HG-U133_Plus_2), ..., BO5649_(HG-U133_Plus_2) [86]
Time period: 2013-03-06 15:28:16 -- 2013-03-06 15:59:18
Total file size: 1111.29MB
RAM: 0.10MB

# Rank-based quantile normalization
qn <- QuantileNormalization(csBC, typesToUpdate="pm",tag="QN");
print(qn);


csN <- process(qn, verbose=verbose);
print(csN);

AffymetrixCelSet:
  Name: responsify
Tags: BC,QN
Path: probeData/responsify,BC,QN/HG-U133_Plus_2
Platform: Affymetrix
Chip type: HG-U133_Plus_2
Number of arrays: 86
Names: BO2098_(HG-U133_Plus_2), BO2196_(HG-U133_Plus_2), BO2229_(HG-U133_Plus_2), ..., BO5649_(HG-U133_Plus_2) [86]
Time period: 2013-03-06 15:28:16 -- 2013-03-06 15:59:18
Total file size: 1111.29MB
RAM: 0.10MB



plm <- RmaPlm(csN)
print(plm)

RmaPlm:
Data set: responsify
Chip type: HG-U133_Plus_2
Input tags: BC,QN
Output tags: BC,QN,RMA
Parameters: {probeModel: chr "pm", shift: num 0, flavor: chr "affyPLM", treatNAsAs: chr "weights"}
Path: plmData/responsify,BC,QN,RMA/HG-U133_Plus_2
RAM: 0.00MB

fit(plm, verbose=verbose) #plm is the probe level estimates


##Quality assessments
#To examine NUSE and RLE plots, do:

qam <- QualityAssessmentModel(plm)
png(file=paste(prefix.out,'Nuse.png', sep = ''), width=1000, height=700)
par(mar=c(15,5,1,1)+0.1)
plotNuse(qam)
dev.off()

png(file=paste(prefix.out,'Rle.png', sep = ''), width=1000, height=700)
par(mar=c(15,5,1,1)+0.1)
plotRle(qam)
dev.off()


cs <- csR
filename <- sprintf("%s,plotDensity.png", getFullName(cs))
png(file=paste(prefix.out,filename, sep = ''), width=640, height=400)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs, lwd=2, ylim=c(0,1))
stext(side=3, pos=0, getFullName(cs))
dev.off()

cs <- csN
filename <- sprintf("%s,plotDensity_after.png", getFullName(cs))
png(file=paste(prefix.out,filename, sep = ''), width=640, height=400)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs, lwd=2, ylim=c(0,0.70))
stext(side=3, pos=0, getFullName(cs))
dev.off()




################################################################
################################################################
#
# Discard bad arrays and extract expression estimates
#
################################################################
################################################################

#Decided to ditch array 26 from NUSE RLE plots

### A couple of MA-plots for further judgements
###

#Get hold of the index vector for the PM probes by:
cs <- csN
cdf <- getCdf(cs);
cells <- getCellIndices(cdf, stratifyBy="pm", unlist=TRUE, useNames=FALSE);
cfR <- getAverageFile(cs); #Average reference "array"


#Keep 19 or 20 (4016bis or 4016)?  
  
cf1 <- getFile(cs, 19);
cf2 <- getFile(cs, 20);
png(file=paste(prefix.out,'MA_4016bis_vs_4016.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf1, cf2, indices=cells);
dev.off()

png(file=paste(prefix.out,'MA_4016bis_vs_pool.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf1, cfR, indices=cells);
dev.off()

png(file=paste(prefix.out,'MA_4016_vs_pool.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf2, cfR, indices=cells);
dev.off()


#Keep 80 or 81 (5416bis or 5416)?  

cf1 <- getFile(cs, 80);
cf2 <- getFile(cs, 81);
png(file=paste(prefix.out,'MA_5419bis_vs_5419.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf1, cf2, indices=cells);
dev.off()

png(file=paste(prefix.out,'MA_5419bis_vs_pool.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf1, cfR, indices=cells);
dev.off()

png(file=paste(prefix.out,'MA_5419_vs_pool.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf2, cfR, indices=cells);
dev.off()


# Compare with a random array

cf1 <- getFile(cs, 1);
png(file=paste(prefix.out,'MA_array1_vs_pool.png', sep = ''), width=640, height=400)
smoothScatterMvsA(cf1, cfR, indices=cells);
dev.off()

#There is possibily some minor difference within the array pairs, but hardly visible.
#We could continue with either, so we choose the bis files.


#To extract the estimates, you can use extractDataFrame() on the ChipEffectSet 
#object that corresponds to the 'plm' object:

ces <- getChipEffectSet(plm)
expr <- extractDataFrame(ces, addNames=TRUE)
expr = expr[,!(names(expr) %in% getNames(ces)[c(20, 26, 81)])]
save(expr = expr, file = paste(prefix.out, 'Expression_estimates.RData', sep = ''))


save.image(paste(prefix.out, 'Expression.RData', sep = '')) 
load(paste(prefix.out, 'Expression.RData', sep = ''))
library(aroma.affymetrix)
date = format(Sys.Date())


dates = getTimestamps(csN)
write.table(dates, paste(prefix.out, 'Array_dates.txt', sep = ''))
################################################################
################################################################
#
# Attach gene names
#
################################################################
################################################################

# Convert the object to a list
xx <- as.list(hgu133plus2ALIAS2PROBE) 
if(length(xx) > 0){
  # Get the probe identifiers for the first two aliases xx[1:2]
  # Get the first one
  xx[[1]]
}


ugp <- getAromaUgpFile(cdf)
data <- readDataFrame(ugp)

library(hgu133plus2.db)
x <- hgu133plus2GENENAME
# Get the probe identifiers that are mapped to a gene name 
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the GENENAME for the first five probes xx[1:5]
  # Get the first one
  xx[[1]]
}

######################################################################################
######################################################################################
END OF FILE








# Shortcut to setup the already background-corrected CEL files
csBC <- AffymetrixCelSet$byName('responsify', tags="BC", cdf=cdf);
print(csBC);
# Shortcut to setup the background-corrected and normalized CEL files
csN <- AffymetrixCelSet$byName('responsify', tags="QN", cdf=cdf);
print(csN);
