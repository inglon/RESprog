################################################################
################################################################
# Script: ProgR_1.1_SNP.R
# Author: Ingrid Lonnstedt
# Date:  15/02/2013
# R version: R 2.15.1
# Details: SNP analysis
# Version 3 has added normal germline SNP arrays.
################################################################
################################################################



################################################################
################################################################
#
# File paths, functions and date
#
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY/AROMA', sep=''))
#setwd(paste(getwd(), '/AROMA', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/reports/SNP_preprocessing/", sep='')

#prefix.prog = paste(getwd(), "/programs/", sep='')
#source(paste(prefix.prog, 'ProgR_0.1_Functions.R', sep=''))
date = format(Sys.Date())
#library(crlmm)
library(aroma.affymetrix)
setOption(aromaSettings, "memory/ram", 200);

#This line assumes the info files have been created (next section):
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))

################################################################
################################################################
#
# Construct info files: Don't run this again
#
################################################################
################################################################

### Array names: This is old and don't ues any more
###
filenames = list.files(prefix.raw)
filenames = filenames[substr(filenames, 1, 1) == 'A']  
filenames = substr(filenames, 1, nchar(filenames, type = 'c')-4)
info = data.frame(arrname = filenames)
info$AROS.ID = substr(info$arrname, 1, 8)
tmp = read.delim(paste(prefix.ann, 'A1761_work sheet_ SNP 6.0_okay_Scan_info.txt', sep = ''),
                 sep = '\t')
info = merge(info, tmp, all = T)[,c('AROS.ID','arrname','Sample.ID')]
info$Sample.ID[info$AROS.ID %in% c('A1974-08','A1974-09','A1974-10')] = 
  c('3860 T','4541 T','4868 T')
#Choose to exclude some arrays after the normalization below:
info$pass = T
info$pass[info$arrname %in% c('A1761-02','A1761-08','A1761-08rh',
                              'A1761-20','A1761-20rh','A1761-26',
                              'A1761-30','A1761-30rh','A1761-50')] = F



#Change names of germline SNP arrays: don't do this again
filenames = list.files(prefix.raw)
filenames = filenames[substr(filenames, 1, 1) == 'A']  
tochange = filenames[substr(filenames, 2, 2) == '2']
for (i in 1:length(tochange)){
  newname = paste(substr(tochange[i], 1, 5), '.', substr(tochange[i], 7, 8), '.CEL', sep = '')
  file.rename(paste(prefix.raw, tochange[i], sep = ''), paste(prefix.raw, newname, sep = ''))
}

### Array names new version
###

info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
  
filenames = list.files(prefix.raw)
filenames = filenames[substr(filenames, 1, 1) == 'A']  
filenames = filenames[substr(filenames, 2, 2) == '2']  
filenames = substr(filenames, 1, nchar(filenames, type = 'c')-4)

moreinfo = data.frame(arrname = filenames, stringsAsFactors = F)
moreinfo$AROS.ID = paste(substr(moreinfo$arrname, 1, 5), '-', substr(moreinfo$arrname, 7, 8),
                         sep = '')

tmp = read.delim(paste(prefix.ann, 'A2163_worksheet_OK_Scaninfo.txt', sep = ''),
                 sep = '\t')
moreinfo = merge(moreinfo, tmp, all = T)[,c('AROS.ID','arrname','Sample.ID')]

info$tissue = 'T'
moreinfo$tissue = 'N'
moreinfo$AROS.ID2 = moreinfo$arrname
moreinfo$arrname2 = moreinfo$AROS.ID
moreinfo$sample = as.integer(substr(moreinfo$Sample.ID, 1, 4))

tmp = subset(moreinfo, select = c('sample', 'AROS.ID'))
names(tmp)[2] = 'myaros'
infotmp = merge(info, tmp, all = T, sort = F)
infotmp$myind = infotmp$ind + 100
infotmp = subset(infotmp, !is.na(myaros), select = c('myaros','myind'))
names(infotmp) = c('AROS.ID','ind')
moreinfo = merge(moreinfo, infotmp, all =T)
moreinfo$sample.name = paste(moreinfo$ind,moreinfo$sample, sep = '.')

info = rbind(info, moreinfo)
info = info[order(info$ind),]
write.table(info, file = paste(prefix.ann, '../../Infokey.txt', sep = ''), 
            quote = F, sep = '\t', row.names = F) 

################################################################
################################################################
#
# Aroma copy numbers init
#
################################################################
################################################################

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
# Don't display too many decimals.'
options(digits=4)

#Make sure annotation files can be found
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

AffymetrixCdfFile:
  Path: annotationData/chipTypes/GenomeWideSNP_6
Filename: GenomeWideSNP_6,Full.cdf
File size: 470.44 MB (493291745 bytes)
Chip type: GenomeWideSNP_6,Full
RAM: 0.00MB
File format: v4 (binary; XDA)
Dimension: 2572x2680
Number of cells: 6892960
Number of units: 1881415
Cells per unit: 3.66
Number of QC units: 4

gi <- getGenomeInformation(cdf)
print(gi)

UgpGenomeInformation:
  Name: GenomeWideSNP_6
Tags: Full,na31,hg19,HB20110328
Full name: GenomeWideSNP_6,Full,na31,hg19,HB20110328
Pathname: annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp
File size: 8.97 MB (9407867 bytes)
RAM: 0.00 MB
Chip type: GenomeWideSNP_6,Full

si <- getSnpInformation(cdf)
print(si)

UflSnpInformation:
  Name: GenomeWideSNP_6
Tags: Full,na31,hg19,HB20110328
Full name: GenomeWideSNP_6,Full,na31,hg19,HB20110328
Pathname: annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl
File size: 7.18 MB (7526452 bytes)
RAM: 0.00 MB
Chip type: GenomeWideSNP_6,Full
Number of enzymes: 2

acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
print(acs)

AromaCellSequenceFile:
  Name: GenomeWideSNP_6
Tags: HB20080710
Full name: GenomeWideSNP_6,HB20080710
Pathname: annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,HB20080710.acs
File size: 170.92 MB (179217531 bytes)
RAM: 0.00 MB
Number of data rows: 6892960
File format: v1
Dimensions: 6892960x26
Column classes: raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw, raw
Number of bytes per column: 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
Footer: <createdOn>20080710 22:47:02 PDT</createdOn><platform>Affymetrix</platform><chipType>GenomeWideSNP_6</chipType><srcFile><filename>GenomeWideSNP_6.probe_tab</filename><filesize>341479928</filesize><checksum>2037c033c09fd8f7c06bd042a77aef15</checksum></srcFile><srcFile2><filename>GenomeWideSNP_6.CN_probe_tab</filename><filesize>96968290</filesize><checksum>3dc2d3178f5eafdbea9c8b6eca88a89c</checksum></srcFile2>
  Chip type: GenomeWideSNP_6
Platform: Affymetrix



################################################################
################################################################
#
# Aroma copy numbers declaring the raw dataset
#
################################################################
################################################################


cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
csR <- AffymetrixCelSet$byName("responsify", cdf=cdf)
print(csR)

AffymetrixCelSet:
  Name: responsify
Tags: 
  Path: rawData/responsify/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2012-10-05 03:43:52 -- 2013-09-14 12:49:35
Total file size: 5666.55MB
RAM: 0.10MB

#Proceed only with arrays that have info$pass == TRUE (skipped this)
#csR = extract(csR, (1:nrow(info))[info$pass])
#print(csR)


#Quality assessment
cs <- csR
filename <- sprintf("%s,plotDensity.png", getFullName(cs))
png(file=paste(prefix.out,filename, sep = ''), width=640, height=400)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs, lwd=2, ylim=c(0,0.40))
stext(side=3, pos=0, getFullName(cs))
dev.off()
#dev.print(png, file=paste(prefix.out,filename, sep = ''), width=640, height=400)
#Figure: The empirical densities for each of the arrays in the data set before any calibration.


ae <- ArrayExplorer(csR) #Only done for 56 arrays
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"))
process(ae, interleaved="auto", verbose=verbose)
print(ae)

ArrayExplorer:
  Name: responsify
Tags: 
  Number of chip types: 1
Number of arrays: 56
Color maps: log2,log2neg,rainbow; log2,log2pos,rainbow
Main path: reports/responsify/raw
RAM: 0.00MB

display(ae)


################################################################
################################################################
#
# Aroma copy numbers step 1:
# Calibration for crosstalk between allele probe pairs
#
################################################################
################################################################

acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
print(acc)

AllelicCrosstalkCalibration:
  Data set: responsify
Input tags: 
  User tags: *
  Asterisk ('*') tags: ACC,ra,-XY
Output tags: ACC,ra,-XY
Number of files: 86 (5666.55MB)
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Algorithm parameters: {rescaleBy: chr "all", targetAvg: num 2200, subsetToAvg: int [1:6584394] 1 2 3 4 5 6 7 8 9 10 ..., mergeShifts: logi TRUE, B: int 1, flavor: chr "sfit", algorithmParameters:List of 3, ..$ alpha: num [1:8] 0.1 0.075 0.05 0.03 0.01 0.0025 0.001 0.0001, ..$ q: num 2, ..$ Q: num 98}
Output path: probeData/responsify,ACC,ra,-XY/GenomeWideSNP_6
Is done: FALSE
RAM: 25.12MB

csC <- process(acc, verbose=verbose)
print(csC)


AffymetrixCelSet:
  Name: responsify
Tags: ACC,ra,-XY
Path: probeData/responsify,ACC,ra,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2013-09-19 20:59:24 -- 2013-09-20 00:18:37
Total file size: 5653.39MB
RAM: 0.10MB

cs <- csC
filename <- sprintf("%s,plotDensity.png", getFullName(cs))
png(paste(prefix.out,filename, sep = ''), width=640, height=400)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs, lwd=2, ylim=c(0,0.40))
stext(side=3, pos=0, getFullName(cs))
dev.off()
#dev.print(png, file=filename, width=640, height=400)
#dev.print(png, file=paste(prefix.out,filename, sep = ''), width=640, height=400)
#Figure: The empirical densities for each of the arrays in the data set after crosstalk calibration.



## Have a look at array 1 allele pairs (tumour sample)
##

filename <- sprintf("%s,plotAllelePairs.png", getFullName(csR))
array <- 1
png(file=paste(prefix.out,filename, array, sep = ''), width=800, height=580)
xlim <- c(-500,15000)
plotAllelePairs(acc, array=array, pairs=1:6, what="input", xlim=xlim/3)
dev.off()
#dev.print(png, file=paste(prefix.out,filename, sep = ''), width=800, height=580)
#Figure: Allele probe pair intensities (PMA,PMB) of array NA06985 for the six nucleotide 
#pairs (A,C), (A,G), (A,T), (C,G), (C,T), and (G,T).  Data shown is before calibration.

filename <- sprintf("%s,plotAllelePairs.png", getFullName(csC))
png(file=paste(prefix.out,filename, sep = ''), width=800, height=580)
plotAllelePairs(acc, array=array, pairs=1:6, what="output", xlim=xlim)
dev.off()
#dev.print(png, file=paste(prefix.out,filename, sep = ''), width=800, height=580)
#Figure: Allele probe pair intensities (PMA,PMB) of array NA06985 for the six nucleotide pairs 
#(A,C), (A,G), (A,T), (C,G), (C,T), and (G,T).  Data shown is after calibration.


## Have a look at array 86 allele pairs (germline sample)
##

filename <- sprintf("%s,plotAllelePairs.png", getFullName(csR))
array <- 86
png(file=paste(prefix.out,filename, array, sep = ''), width=800, height=580)
xlim <- c(-500,15000)
plotAllelePairs(acc, array=array, pairs=1:6, what="input", xlim=xlim/3)
dev.off()
#dev.print(png, file=paste(prefix.out,filename, sep = ''), width=800, height=580)
#Figure: Allele probe pair intensities (PMA,PMB) of array NA06985 for the six nucleotide 
#pairs (A,C), (A,G), (A,T), (C,G), (C,T), and (G,T).  Data shown is before calibration.

filename <- sprintf("%s,plotAllelePairs.png", getFullName(csC))
png(file=paste(prefix.out,filename, array, sep = ''), width=800, height=580)
plotAllelePairs(acc, array=array, pairs=1:6, what="output", xlim=xlim)
dev.off()




## Allele pair figures for all arrays (only done for 62 tumour arrays)
##

xlim <- c(-500,15000)
pdf(paste(prefix.out, 'allAllelePairs.pdf', sep=''), paper = 'a4') 
for (array in 1:62){
  plotAllelePairs(acc, array=array, pairs=1:6, what="input", xlim=xlim/3)
  plotAllelePairs(acc, array=array, pairs=1:6, what="output", xlim=xlim)  
}
dev.off()

################################################################
################################################################
#
# Aroma copy numbers step 2:
# Normalization for nucleotide-position probe sequence effects
#
################################################################
################################################################

bpn <- BasePositionNormalization(csC, target="zero")
#By using argument target="zero", no reference is required.  
#Otherwise, the average file will be used as the reference.
print(bpn)


print(csN)
BasePositionNormalization:
  Data set: responsify
Input tags: ACC,ra,-XY
User tags: *
  Asterisk ('*') tags: BPN,-XY
Output tags: ACC,ra,-XY,BPN,-XY
Number of files: 86 (5653.39MB)
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Algorithm parameters: {unitsToFit: int [1:1784456] 1 2 3 4 5 6 7 8 9 10 ..., typesToFit: chr "pm", unitsToUpdate: int [1:1881415] 1 2 3 4 5 6 7 8 9 10 ..., typesToUpdate: chr "pm", shift: num 0, cellsToFit: int [1:6527163] 3466655 3463975 3461295 3458615 3455935 3453255 3450575 3447895 3445215 3442535 ..., cellsToUpdate: int [1:6835685] 3466655 3463975 3461295 3458615 3455935 3453255 3450575 3447895 3445215 3442535 ..., target: chr "zero", model: chr "smooth.spline", df: int 5}
Output path: probeData/responsify,ACC,ra,-XY,BPN,-XY/GenomeWideSNP_6
Is done: FALSE
RAM: 50.98MB

csN <- process(bpn, verbose=verbose) # In transform(y) : NaNs produced
print(csN)


AffymetrixCelSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY
Path: probeData/responsify,ACC,ra,-XY,BPN,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2013-09-19 20:59:24 -- 2013-09-20 00:18:37
Total file size: 5653.39MB
RAM: 0.10MB



cs <- csN
filename <- filename <- sprintf("%s,plotDensity.png", getFullName(cs))
png(file=paste(prefix.out,filename, sep = ''), width=640, height=400)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs, lwd=2, ylim=c(0,0.40))
stext(side=3, pos=0, getFullName(cs))
dev.off()
# dev.print(png, file=filename, width=640, height=400)
#Figure: The empirical densities for each of the arrays in the data set after 
#crosstalk calibration and nucleotide-position normalization.


## Look at array 1 allele pairs
##
array <- 1
xlim <- c(-500,15000)
filename <- sprintf("%s,plotAllelePairs", array,".png", getFullName(csN))
acc2 <- AllelicCrosstalkCalibration(csN)
png(file=paste(prefix.out,filename, sep = ''), width=800, height=580)
plotAllelePairs(acc2, array=array, pairs=1:6, what="input", xlim=1.5*xlim)
dev.off()
#dev.print(png, file=filename, width=800, height=580)


## Look at array 1 allele pairs
##
array <- 86
xlim <- c(-500,15000)
filename <- sprintf("%s,plotAllelePairs", array,".png", getFullName(csN))
acc2 <- AllelicCrosstalkCalibration(csN)
png(file=paste(prefix.out,filename, sep = ''), width=800, height=580)
plotAllelePairs(acc2, array=array, pairs=1:6, what="input", xlim=1.5*xlim)
dev.off()

#Figure: Allele probe pair intensities (PMA,PMB) of array NA06985 for the six nucleotide 
#pairs (A,C), (A,G), (A,T), (C,G), (C,T), and (G,T).  Data shown is after crosstalk calibration
#and nucleotide-position normalization.  Note how the heterozygote arms are along the diagonals, 
#that is, there is a balance in the allele A and allele B signal for heterozygotes.  
#This is (on purpose) not corrected for in the allelic crosstalk calibration.

## Allele pair figures for all arrays (only done with tumour arrays)
##

array <- 1
xlim <- c(-500,15000)
acc2 <- AllelicCrosstalkCalibration(csN)
pdf(paste(prefix.out, getFullName(csN), ',allAllelePairs.pdf', sep=''), paper = 'a4') 
for (array in 1:62){
  plotAllelePairs(acc2, array=array, pairs=1:6, what="input", xlim=1.5*xlim)
}
dev.off()


################################################################
################################################################
#
# Aroma copy numbers step 3:
# Probe summarization
#
################################################################
################################################################

#Next we summarize the probe level data unit by unit.  For SNPs we have the option 
#to model either the total CN signals (combineAlleles=TRUE) or 
#allele-specific signals (combineAlleles=FALSE).  Here we fit allele specific CN signals.
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE)
print(plm)

RmaCnPlm:
  Data set: responsify
Chip type: GenomeWideSNP_6,Full
Input tags: ACC,ra,-XY,BPN,-XY
Output tags: ACC,ra,-XY,BPN,-XY,RMA
Parameters: {probeModel: chr "pm", shift: num 0, flavor: chr "affyPLM", treatNAsAs: chr "weights", mergeStrands: logi TRUE, combineAlleles: logi FALSE}
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA/GenomeWideSNP_6
RAM: 0.00MB

if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=verbose)
  str(units)
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...
  
  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=verbose)
  str(units)
}

ces <- getChipEffectSet(plm)
print(ces)
#Pasted to here 20 sept 2013

CnChipEffectSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2013-09-20 10:43:56 -- 2013-09-20 10:44:33
Total file size: 2317.63MB
RAM: 0.14MB
Parameters: {}





################################################################
################################################################
#
# Aroma copy numbers step 4:
# Normalization for PCR fragment-length effects
#
################################################################
################################################################

#Similarly to how we normalized for the probe-sequence effects, 
#we will here normalize for PCR fragment-length effects by using a "zero" 
#target.  This will avoid using the average (chip effects) as a reference.  
#Thus, this step is also truly single-array by nature.

fln <- FragmentLengthNormalization(ces, target="zero")
print(fln)

FragmentLengthNormalization:
  Data set: responsify
Input tags: ACC,ra,-XY,BPN,-XY,RMA
User tags: *
  Asterisk ('*') tags: FLN,-XY
Output tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Number of files: 86 (2317.63MB)
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Algorithm parameters: {subsetToFit: chr "-XY", lengthRange: NULL, onMissing: chr "median", .target: chr "zero", shift: num 0}
Output path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Is done: FALSE
RAM: 13.62MB

cesN <- process(fln, verbose=verbose)
print(cesN)

CnChipEffectSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2013-09-20 10:43:56 -- 2013-09-20 10:44:33
Total file size: 2317.63MB
RAM: 0.14MB
Parameters: {}








################################################################
################################################################
#
# Quality assessments
# 
# 
################################################################
################################################################

ae <- ArrayExplorer(csR) #
setColorMaps(ae, "sqrt,yellow")
process(ae, verbose=TRUE)
print(ae)

ArrayExplorer:
  Name: responsify
Tags:
  Number of chip types: 1
Number of arrays: 86
Color maps: sqrt,yellow
Main path: reports/responsify/raw
RAM: 0.00MB

#Before normalization for PCR fragment-length effects
qam <- QualityAssessmentModel(plm)
print(qam)

QualityAssessmentModel:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,QC
Path: qcData/responsify,ACC,ra,-XY,BPN,-XY,RMA,QC/GenomeWideSNP_6
Chip-effect set:
  CnChipEffectSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Number of arrays: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Time period: 2013-09-20 10:43:56 -- 2013-09-20 10:44:33
Total file size: 2317.63MB
RAM: 0.14MB
Parameters: {}
RAM: 0.00MB

png(file=paste(prefix.out,'Nuse.png', sep = ''), width=1000, height=700)
par(mar=c(8,5,1,1)+0.1)
plotNuse(qam)
dev.off()

png(file=paste(prefix.out,'Rle.png', sep = ''), width=1000, height=700)
par(mar=c(8,5,1,1)+0.1)
plotRle(qam)
dev.off()


#Looked at array dates, batches can be found from array names.
dates = getTimestamps(csN)
#Pasted to here friday but have not pasted all the printouts.


################################################################
################################################################
#
# Aroma copy numbers step 5:
# CalMaTe post-normalization with library calmate
#
################################################################
################################################################
library("calmate");
dsNList <- exportTotalAndFracB(cesN, verbose=verbose);

names <- getNames(dsNList$total);
idxsN <- which(substr(names, 2, 2) == '2' & names != 'A2163.02');


#Next, when setting up the CalMaTe model, we specify the reference samples using argument 'references' as follows:
cmtN <- CalMaTeCalibration(dsNList, tags=c("*", "normalReferences"), references=idxsN);
print(cmtN);

CalMaTeCalibration:
  Data sets (2):
  <Total>:
  AromaUnitTotalCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Number of files: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Total file size: 617.26 MB
RAM: 0.11MB
<FracB>:
  AromaUnitFracBCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Number of files: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Total file size: 617.26 MB
RAM: 0.10MB
Number of arrays: 86
Number of references: 29 (33.72%)
Additional parameters: [1] {references: int [1:29] 57 59 60 61 62 63 64 65 66 67 ..., flavor: chr "v2"}

dsCList <- process(cmtN, verbose=verbose);
print(dsCList);

$total
AromaUnitTotalCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences
Number of files: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences/GenomeWideSNP_6
Total file size: 617.26 MB
RAM: 0.11MB

$fracB
AromaUnitFracBCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences
Number of files: 86
Names: A1761.01, A1761.02rh, A1761.03, ..., A2163.30 [86]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2,normalReferences/GenomeWideSNP_6
Total file size: 617.26 MB
RAM: 0.10MB

#List of sample which look like we should be able to get CNs:
out = subset(info, select = c('ind','AROS.ID','sample'))
out$analyzable = 'Unlikely'
out$analyzable[c(1:2, 5, 9, 11:13, 15, 18:21, 23:27, 29:35, 38:40, 44:47, 49, 51:52,
                 55:56)] = 'OK'
out$analyzable[c(3, 4, 8, 14, 16, 28, 37, 41, 43, 50, 53, 54)] = 'Possibly'
normals = out[57:86,]
normals = normals[normals$ind != 154,]
myind = normals$ind-100
out = out[1:56,]
out$hasNormal = F
out$hasNormal[out$ind %in% myind]= T
out$analyzableAndNormal = F
out$analyzableAndNormal[out$analyzable == 'OK' & out$hasNormal == T] = T


################################################################
################################################################
#
# Segmentation with library PSCBS
# 
# 
################################################################
################################################################
prefix.out =  "/wehisan/home/allstaff/l/lonnstedt/RESPONSIFY/AROMA/reports/SNP_segmentation/"
library(PSCBS)
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

#extract (total, beta) for samples
dataC = extractPSCNArray(dsCList$total)

#Look at betaT
pdf(paste(prefix.out, 'BetaTwithNormals.pdf', sep=''), paper = 'a4') 
par(mfrow = c(2,1))
for (arr in 1:dim(dataC)[3]){
  
  #Prepare data for this sample
  samp = as.character(info$arrname[arr])
  betaT = dataC[,'fracB',samp]
  betaT = betaT[!is.na(betaT)]
  lims = quantile(betaT, na.rm = T, probs = c(.025, .975))
  toplot = betaT[betaT>lims[1] & betaT<lims[2]]
  plotDensity(toplot, xlim = c(0, 1), main = arr)
  grid()
}  
dev.off()

tauAs = numeric(56)
tauAs[c(5, 12, 14:15, 19:23, 27:28, 32:34, 37:38, 40:42, 56)] = .2
tauAs[c(1:4, 6,11, 13, 16:17, 24, 30, 36, 46)] = .25
tauAs[c(7:8, 29, 31, 35, 39, 47:48, 51)] = .225
tauAs[c( 9:10, 18, 43:45, 55)] = .15
tauAs[c(25 )] = .175
tauAs[c( 26)] = .3

normix = which(info$tissue == 'N' & info$AROS.ID2 != 'A2163.02')
tumix = which(info$tissue == 'T')
myNList <- lapply(dsCList, extract, normix);
yR <- getAverageFile(myNList$total, verbose = T);

### Segmentation for samples with no normal
###

for (arr in 1:56){ #This should have been done only for samples without normals
  
  #Prepare data for this sample
  samp = as.character(info$arrname[arr])
  fitname = paste('fit', arr, sep = '')
  yT = dataC[,'total',samp]
  if ((arr+100) %in% info$ind & arr != 54){
    sampN = as.character(info$arrname[info$ind == (arr + 100)])
    yN = dataC[,'total',sampN]
    CT = yT/yN
  } else   CT = (yT/yR[,1])[,1]
  betaT = dataC[,'fracB',samp]
  ugp = getAromaUgpFile(dsNList$total)
  chromosome = ugp[,1,drop = T]
  x = ugp[, 2, drop = T]
  df = data.frame(chromosome = chromosome, x=x, CT=CT, betaT = betaT)
  df = dropSegmentationOutliers(df)
  gaps = findLargeGaps(df, minLength = 1e+06)
  knownSegments = gapsToSegments(gaps)
  
  ## 'betaT' is already normalized using TumorBoost => tbn=FALSE.
  fit <- segmentByNonPairedPSCBS(df, knownSegments=knownSegments, tauA = tauAs[arr],
                             avgDH="median", seed=0xBEEF, verbose=-10);  
  
  ## Estimate global background level in [0,1] (due to normal
  #contamination and more)
  kappa <- try(estimateKappa(fit, verbose=-10), silent = T)
  
  ## Call allelic balance
  ## (a) Estimate DH threshold for calling AB
  if (substr(kappa[1], 1, 5) != 'Error'){
    deltaAB <- estimateDeltaAB(fit, kappa=kappa); # If skipped, will be done internally
    
    ## (b) Call AB based on bootstrapped segment DH levels
    fit <- callAB(fit, delta=deltaAB, verbose=-100);
  } else fit = callAB(fit, verbose = -100)
  
  ## Call loss of heterozygosity (LOH)
  ## (a) Estimate C_1 threshold for calling LOH
  deltaLOH <- estimateDeltaLOH(fit); # If skipped, will be done internally
  
  ## (b) Call LOH based on bootstrapped segment C_1 levels
  fit <- callLOH(fit, delta=deltaLOH, verbose=-10);
  
  ## Call NTCN
  ## (a) Estimate the threshold for calling neutral TCN segments
  ##     By shrinking 'scale', more segments will be non-NTCN.
  if (substr(kappa[1], 1, 5) != 'Error'){
    deltaCN <- estimateDeltaCN(fit, scale=1.0, kappa=kappa);
    fit <- callNTCN(fit, delta=deltaCN, verbose=-20);
    
    ## (b) Call NTCN based on bootstrapped segment TCN levels
  } else {
    deltaCN = try(estimateDeltaCN(fit, scale = 1.0), silent = T)
    if (substr(deltaCN[1], 1, 5) != 'Error') fit <- callNTCN(fit, delta=deltaCN, verbose=-20) else {
      tmp = try(callNTCN(fit, verbose = -20), silent = T)
      if (substr(tmp[1], 1, 5) != 'Error') fit = tmp
    }
    
  }
  
  #  toPNG(fitname, tags="tracks,avgDH=median,AB+LOH+NTCN", aspectRatio=0.6, {
  #    plotTracks(fit, chromosomes=1:23);
  #  });
  
  filename <- paste(prefix.out, 'CBSfitsWithNormals/', fitname, '.RData', sep = '');
  saveObject(fit, file=filename);
  png(file=paste(prefix.out, 'CBSfitsWithNormals/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
  rm(list = c('fit','kappa','deltaAB','deltaLOH','deltaCN'))
}  

save.image(paste(prefix.out, 'SNP.RData', sep = '')) 


### Segmentation for samples with normals matched pair
###

for (arr in setdiff(tumix, c(21))){ 
  
  if ((arr+100) %in% info$ind & arr != 54){
    
    #Prepare data for this sample
    samp = as.character(info$arrname[arr])
    fitname = paste('fit', arr, sep = '')
    yT = dataC[,'total',samp]
    sampN = as.character(info$arrname[info$ind == (arr + 100)])
    yN = dataC[,'total',sampN]
    CT = yT/yN
    betaT = dataC[,'fracB',samp]
    betaN = dataC[,'fracB',sampN]  
    ugp = getAromaUgpFile(dsNList$total)
    chromosome = ugp[,1,drop = T]
    x = ugp[, 2, drop = T]
    df = data.frame(chromosome = chromosome, x=x, CT=CT, betaT = betaT, betaN = betaN)
    df = dropSegmentationOutliers(df)
    gaps = findLargeGaps(df, minLength = 1e+06)
    knownSegments = gapsToSegments(gaps)
    
    ## 'betaT' is already normalized using TumorBoost => tbn=FALSE.
    fit <- segmentByPairedPSCBS(df, knownSegments=knownSegments, tauA = tauAs[arr],
                                avgDH="median", seed=0xBEEF, verbose=-10);  
    
    try({
    ## Estimate global background level in [0,1] (due to normal
    #contamination and more)
    kappa <- try(estimateKappa(fit, verbose=-10), silent = T)
    
    ## Call allelic balance
    ## (a) Estimate DH threshold for calling AB
    if (substr(kappa[1], 1, 5) != 'Error'){
      deltaAB <- estimateDeltaAB(fit, kappa=kappa); # If skipped, will be done internally
      
      ## (b) Call AB based on bootstrapped segment DH levels
      fit <- callAB(fit, delta=deltaAB, verbose=-100);
    } else fit = callAB(fit, verbose = -100)
    
    ## Call loss of heterozygosity (LOH)
    ## (a) Estimate C_1 threshold for calling LOH
    deltaLOH <- estimateDeltaLOH(fit); # If skipped, will be done internally
    
    ## (b) Call LOH based on bootstrapped segment C_1 levels
    fit <- callLOH(fit, delta=deltaLOH, verbose=-10);
    
    ## Call NTCN
    ## (a) Estimate the threshold for calling neutral TCN segments
    ##     By shrinking 'scale', more segments will be non-NTCN.
    if (substr(kappa[1], 1, 5) != 'Error'){
      deltaCN <- estimateDeltaCN(fit, scale=1.0, kappa=kappa);
      fit <- callNTCN(fit, delta=deltaCN, verbose=-20);
      
      ## (b) Call NTCN based on bootstrapped segment TCN levels
    } else {
      deltaCN = try(estimateDeltaCN(fit, scale = 1.0), silent = T)
      if (substr(deltaCN[1], 1, 5) != 'Error') fit <- callNTCN(fit, delta=deltaCN, verbose=-20) else {
        tmp = try(callNTCN(fit, verbose = -20), silent = T)
        if (substr(tmp[1], 1, 5) != 'Error') fit = tmp
      }
      
    }
    }, silent = T)
    
    filename <- paste(prefix.out, 'CBSfitsWithNormals/', fitname, '.RData', sep = '');
    saveObject(fit, file=filename);
    png(file=paste(prefix.out, 'CBSfitsWithNormals/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
    plotTracks(fit, chromosomes=1:23)
    title(main = samp, outer = T)
    dev.off()   
    rm(list = c('fit','kappa','deltaAB','deltaLOH','deltaCN'))
  }
}  

save.image(paste(prefix.out, 'SNP.RData', sep = '')) 

#Pasted to here


### Print segmentation figures without NTCN
###
load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

##NOTE: NTCN calls are not optimal, so do not use them!
##Print segmentations without the NTCN cutoffs:
for (arr in tumix){
  
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfitsWithNormals/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  fit$output$ntcnCall = NULL
  png(file=paste(prefix.out, 'CBSfitsWithNormals/CNsegments_samp', arr, 'noNTCN.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
} 


###Number of segments in each sample
###
nsegs = NULL
for (arr in tumix){
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfitsWithNormals/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  nsegs = c(nsegs, nrow(fit$output))
} 
names(nsegs) = 1:56

nsegs

#This was with all segmentations based on nonPaired function:
1     2     3     4     5     6     7     8     9    10    11    12    13 
664   808   493  1003   819   516   835   616   669   442   374   646   862 
14    15    16    17    18    19    20    21    22    23    24    25    26 
374   939   601   428   774   322   278   918   384 17735   914   811   732 
27    28    29    30    31    32    33    34    35    36    37    38    39 
609   757   320   397   321   668   935   596   562   584   811   727  7411 
40    41    42    43    44    45    46    47    48    49    50    51    52 
512   507   161   892   303   779   523   818   569   714   540   735   686 
53    54    55    56 
518   395   213   483 

#These are final:
1     2     3     4     5     6     7     8     9    10    11    12    13 
664   808   493  1008   819   516   835   616   684   451   393   646   862 
14    15    16    17    18    19    20    21    22    23    24    25    26 
380   939   604   437   774   322   279   925   380 17735   905   815   751 
27    28    29    30    31    32    33    34    35    36    37    38    39 
609   757   316   398   322   673   966   596   562   584   811   727  7426 
40    41    42    43    44    45    46    47    48    49    50    51    52 
512   518  1691   901   365   779   519   818   569   728   541   735   696 
53    54    55    56 
518   395  1962   490 

################################################################
################################################################
#
# Data for TAPS
# 
# 
################################################################
################################################################
setwd(paste(getwd(), '/AROMA', sep='')) 
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.out = paste(getwd(), "/reports/SNP_segmentation/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.taps = paste(getwd(), '../TAPS/TAPSWithNormals/', sep = '/')

load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

### Array names
###
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))

# Identification of units in Chr 2:81-86MB and their positions
cdf <- getCdf(cesN) 
gi <- getGenomeInformation(cdf)
pos <- getAromaUgpFile(gi)


for (arr in setdiff(1:56, c(23, 39))){
  #Prepare data for this sample
  samp = as.character((info$arrname)[arr])
  
  yT = dataC[,'total',samp]
  if ((arr+100) %in% info$ind & arr != 54){
    sampN = as.character(info$arrname[info$ind == (arr + 100)])
    yN = dataC[,'total',sampN]
    CT = yT/yN
  } else   CT = (yT/yR[,1])[,1]
  
  log2 = log2( CT) #Warnings OK
  betaT = dataC[,'fracB',samp]
  probes = data.frame(Chromosome = pos[,1], Start = pos[,2], Value = log2)
  colnames(probes) = c('Chromosome','Start','Value')
  probes = subset(probes, !is.na(probes$Chromosome) & !is.na(probes$Value) & 
                    !is.na(probes$Start) & !(probes$Value %in% c(Inf, -Inf)))
  probes$Chromosome = paste("chr", probes$Chromosome, sep = '')
  snps = data.frame(Chromosome = pos[,1], Start = pos[,2], Value = betaT)
  colnames(snps) = c('Chromosome','Start','Value')
  snps = subset(snps, !is.na(snps$Chromosome) & !is.na(snps$Value) & !is.na(snps$Start) & 
                  !(snps$Value %in% c(Inf, -Inf)))
  snps$Chromosome = paste("chr", snps$Chromosome, sep = '')
  probes$Chromosome = as.character(probes$Chromosome)
  probes$Chromosome[probes$Chromosome == "chr23"] = "chrX"
  probes$Chromosome[probes$Chromosome == "chr24"] = "chrY"
  snps$Chromosome = as.character(snps$Chromosome)
  snps$Chromosome[snps$Chromosome == "chr23"] = "chrX"
  snps$Chromosome[snps$Chromosome == "chr24"] = "chrY"
  probes = subset(probes, Chromosome != "chr25")
  snps = subset(snps, Chromosome != "chr25")
  
  sample.name <- paste(arr, (substr(as.character(info$Sample.ID), 1, 4))[arr], sep = '.')
  mydir = paste(prefix.taps, sample.name, sep ='')
  dir.create(mydir, showWarnings = F, recursive = T)
  write.table(probes, paste(mydir, '/probes.txt', sep = ''), row.names=F, col.names=T, quote = F,
              sep = '\t')
  write.table(snps, paste(mydir, '/snps.txt', sep = ''), row.names=F, col.names=T, quote = F,
              sep = '\t')
  
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, '../SNP_segmentation/CBSfitsWithNormals/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  segments = subset(fit$output, !is.na(chromosome) & !is.na(tcnMean) & !is.na(tcnStart) & !is.na(tcnEnd))
  segments$chromosome = paste("chr", segments$chromosome, sep = '')
  segments = data.frame(Chromosome = segments$chromosome, Start = segments$tcnStart,
                        End = segments$tcnEnd, Value = log2(segments$tcnMean/2))
  segments$Chromosome = as.character(segments$Chromosome)
  segments$Chromosome[segments$Chromosome == "chr23"] = "chrX"
  segments$Chromosome[segments$Chromosome == "chr24"] = "chrY"
  segments = subset(segments, !(segments$Value %in% c(Inf, -Inf)))
  segments = subset(segments, Chromosome != "chr25")
  write.table(segments, paste(mydir, '/segments.txt', sep = ''), row.names=F, col.names=T, quote = F,
              sep = '\t')
}


################################################################
################################################################
#
# Integer copy numbers from TAPS
# 
# 
################################################################
################################################################
Workflow:
1. From the folder containing your samples (sample folders) run TAPS_plot().
2. Investigate the scatter plots generated in your sample folders.
3. To proceed with copy number calls, find and open the file "SampleData.txt".
4. For each sample, enter an interpretation of Log-Ratio @ copy number 2 ("cn2"), the difference in Log-Ratio to a deletion ("delta") and the allelic imbalance ratio of CNNLOH ("loh"). Save the file.
5. Run TAPS_call().
6. Inspect the karyotype_check images, and the new chromosome-wise images.
7. If all looks reasonable, you will find good copy number estimates in 'Copynumbers.csv'.
9. Be wary of the result on sex chromosomes which may be difficult to auto-interpret.
10. Watch all images for signs of segmentation failure and tumor cell heterogeneity.

library(TAPS)
prefix.taps = paste(getwd(), '/TAPS/TAPSWithNormals/', sep = '/')
TAPS_plot(directory = prefix.taps)

#To do TAPS_call:
library(TAPS)

prefix.taps = paste(getwd(), 'TAPS/', sep = '/')
myload.txt <- function(file, ...) { 
  read.csv(file,sep='\t', ...) 
}
deogram = TAPS:::getIdeogram

myTAPS_call <- function(directory=NULL,#xlim=c(-1,1),ylim=c(0,1),
                      minseg=1,maxCn=12) {
  ## TAPS_call outputs the total and minor allele copy numbers of all segments as a text file, and as images for visual confirmation.
  ## sampleInfo_TAPS.txt must be present in each sample folder. If TAPS_plot could not make a good guess of the Log-R of copy number 2 
  ## and the Log-R difference to a deletion, you must interpret the scatter plots and edit sampleInfo_TAPS.txt.
  if (is.null(directory))
  {
    cat("No directory supplied, using working directory.")
    directory = "."
    #cat("You have not assigned a directory containing one or more folders of samples for TAPS_call to execute. \n")
    #cat("Example: \"/user/mysamples/\" or, to run it in your current working directory, \".\" \n")
    #directory = readline("Please supply such a directory now: ")
  }
  
  
  setwd(directory)
  #subs <- getSubdirs()
  
  if (length(grep('SampleData.txt',dir()))==1)
  {
    sampleData <- myload.txt('SampleData.txt', colClasses = c('character', rep('numeric', 5), 'character'))
  }
  else
  {
    sampleData <- myload.txt('../SampleData.txt', colClasses = c('character', rep('numeric', 5), 'character'))
  }
  subs=as.character(sampleData$Sample)
  
  if (is.null(subs)) {
    subs=thisSubdir()
    setwd('..')
  }
  for (i in 1:length(subs)) if (sampleData$calculate.copynumbers[i]=='yes') {
    setwd(subs[i])
    name <- subs[i]
    sampleInfo <- sampleData[sampleData$Sample==subs[i],]
    if (nrow(sampleInfo)==1) {
      
      cat(' ..loading', subs[i])
      Log2 <- TAPS:::readLog2()
      alf <- TAPS:::readAlf(localDir)
      segments <- TAPS:::readSegments()
      
      #Some samples throw NA values, we simply remove these.
      Log2=Log2[!is.nan(Log2$Value),]
      Log2=Log2[!is.na(Log2$Value),]
      
      alf=alf[!is.nan(alf$Value),]
      alf=alf[!is.na(alf$Value),]
      
      segments <- segments[!is.nan(segments$Value),]
      segments <- segments[!is.na(segments$Value),]    
      
      segments$Value <- segments$Value-median(Log2$Value) 
      Log2$Value <- Log2$Value-median(Log2$Value)
      
      cat(' ..processing.\n')
      
      load('allRegions.Rdata')                            ## These were prepared in TAPS_plot
      #allRegions <- makeRegions(Log2, alf, segments)
      
      ## estimates the Log-R and Allelic Imbalance Ration of all variants up to maxCn
      t <- TAPS:::findCNs(Log2,alf,allRegions,dmin=0.9,maxCn=maxCn,ceiling=1,sampleInfo=sampleInfo) 
      
      u <- TAPS:::setCNs(allRegions,t$int,t$ai,maxCn)            ## Assigns copy number variant for all segments
      allRegions$regions <- u$regions
      ## adjacent segments with idendical copy number are merged (except over centromere) and all are saved to a text file
      TAPS:::save.txt(u$regions,file=paste(name,'_segmentCN.txt',sep='')) 
      regions=allRegions$regions
      save(t,regions,file="regions_t.Rdata")
      
      TAPS:::karyotype_check(regions$Chromosome,regions$Start,regions$End,regions$log2,regions$imba,regions$Cn,regions$mCn,t,ideogram=NULL,name=name)
      
      TAPS:::karyotype_chromsCN(regions$Chromosome,regions$Start,regions$End,regions$log2,
                         regions$imba,regions$Cn,regions$mCn,ideogram=NULL,
                         as.character(Log2$Chromosome),Log2$Start,Log2$Value,as.character(alf$Chromosome),
                         alf$Start,alf$Value,t,name=name,xlim=c(-1,1),ylim=c(0,1))
      
      cat('..done\n')
      sampleData$calculate.copynumbers[i] <- 'done'
    } else cat('Skipped',name,'\n')
    
    setwd('..')
  }
  TAPS:::save.txt(sampleData,file='SampleData.txt')
}

myTAPS_call(directory = prefix.taps) ##We have not done this!

################################################################
################################################################
