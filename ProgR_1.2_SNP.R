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

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/output/SNP/", sep='')

#prefix.prog = paste(getwd(), "/programs/", sep='')
#source(paste(prefix.prog, 'ProgR_0.1_Functions.R', sep=''))
date = format(Sys.Date())
#library(crlmm)
library(aroma.affymetrix)
setOption(aromaSettings, "memory/ram", 200);

### Array names
###
filenames = list.files(prefix.raw)
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
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2012-10-05 03:11:05 -- 2013-04-14 15:44:08
Total file size: 4283.17MB
RAM: 0.08MB

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
Number of files: 65 (4283.17MB)
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Algorithm parameters: {rescaleBy: chr "all", targetAvg: num 2200, subsetToAvg: int [1:6584394] 1 2 3 4 5 6 7 8 9 10 ..., mergeShifts: logi TRUE, B: int 1, flavor: chr "sfit", algorithmParameters:List of 3, ..$ alpha: num [1:8] 0.1 0.075 0.05 0.03 0.01 0.0025 0.001 0.0001, ..$ q: num 2, ..$ Q: num 98}
Output path: probeData/responsify,ACC,ra,-XY/GenomeWideSNP_6
Is done: TRUE
RAM: 51.21MB

csC <- process(acc, verbose=verbose)
print(csC)


AffymetrixCelSet:
  Name: responsify
Tags: ACC,ra,-XY
Path: probeData/responsify,ACC,ra,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2013-03-06 11:54:45 -- 2013-05-01 12:34:19
Total file size: 4272.91MB
RAM: 0.08MB

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



## Have a look at array 1 allele pairs
##

filename <- sprintf("%s,plotAllelePairs.png", getFullName(csR))
array <- 1
png(file=paste(prefix.out,filename, sep = ''), width=800, height=580)
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

## Allele pair figures for all arrays (only done for 62 arrays)
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
Number of files: 65 (4272.91MB)
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Algorithm parameters: {unitsToFit: int [1:1784456] 1 2 3 4 5 6 7 8 9 10 ..., typesToFit: chr "pm", unitsToUpdate: int [1:1881415] 1 2 3 4 5 6 7 8 9 10 ..., typesToUpdate: chr "pm", shift: num 0, cellsToFit: int [1:6527163] 3466655 3463975 3461295 3458615 3455935 3453255 3450575 3447895 3445215 3442535 ..., cellsToUpdate: int [1:6835685] 3466655 3463975 3461295 3458615 3455935 3453255 3450575 3447895 3445215 3442535 ..., target: chr "zero", model: chr "smooth.spline", df: int 5}
Output path: probeData/responsify,ACC,ra,-XY,BPN,-XY/GenomeWideSNP_6
Is done: TRUE
RAM: 50.98MB

csN <- process(bpn, verbose=verbose) # In transform(y) : NaNs produced
print(csN)

AffymetrixCelSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY
Path: probeData/responsify,ACC,ra,-XY,BPN,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2013-03-06 11:54:45 -- 2013-05-01 12:34:19
Total file size: 4272.91MB
RAM: 0.08MB



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
filename <- sprintf("%s,plotAllelePairs.png", getFullName(csN))
acc2 <- AllelicCrosstalkCalibration(csN)
png(file=paste(prefix.out,filename, sep = ''), width=800, height=580)
plotAllelePairs(acc2, array=array, pairs=1:6, what="input", xlim=1.5*xlim)
dev.off()
#dev.print(png, file=filename, width=800, height=580)

#Figure: Allele probe pair intensities (PMA,PMB) of array NA06985 for the six nucleotide 
#pairs (A,C), (A,G), (A,T), (C,G), (C,T), and (G,T).  Data shown is after crosstalk calibration
#and nucleotide-position normalization.  Note how the heterozygote arms are along the diagonals, 
#that is, there is a balance in the allele A and allele B signal for heterozygotes.  
#This is (on purpose) not corrected for in the allelic crosstalk calibration.

## Allele pair figures for all arrays
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
RAM: 0.01MB

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

CnChipEffectSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2013-03-25 09:31:57 -- 2013-05-01 15:10:16
Total file size: 1751.70MB
RAM: 0.11MB
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
Number of files: 65 (1751.70MB)
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
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2013-03-25 09:31:57 -- 2013-05-01 15:10:16
Total file size: 1751.70MB
RAM: 0.11MB
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
Number of arrays: 56
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
Number of arrays: 65
Names: A1761-01, A1761-02, A1761-02rh, ..., A1974-10 [65]
Time period: 2013-03-25 09:31:57 -- 2013-05-01 15:10:16
Total file size: 1751.70MB
RAM: 0.11MB
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

#Choose to exclude some arrays:
info$pass = T
info$pass[info$arrname %in% c('A1761-02','A1761-08','A1761-08rh',
                              'A1761-20','A1761-20rh','A1761-26',
                              'A1761-30','A1761-30rh','A1761-50')] = F
cesN = extract(cesN, (1:nrow(info))[info$pass])
print(cesN)

CnChipEffectSet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Path: plmData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Platform: Affymetrix
Chip type: GenomeWideSNP_6,Full,monocell
Number of arrays: 56
Names: A1761-01, A1761-02rh, A1761-03, ..., A1974-10 [56]
Time period: 2013-03-25 09:31:57 -- 2013-05-01 15:10:16
Total file size: 1509.16MB
RAM: 0.09MB
Parameters: {}

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
cmt <- CalMaTeCalibration(dsNList);
print(cmt);


CalMaTeCalibration:
  Data sets (2):
  <Total>:
  AromaUnitTotalCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Number of files: 56
Names: A1761-01, A1761-02rh, A1761-03, ..., A1974-10 [56]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Total file size: 401.94 MB
RAM: 0.07MB
<FracB>:
  AromaUnitFracBCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY
Number of files: 56
Names: A1761-01, A1761-02rh, A1761-03, ..., A1974-10 [56]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY/GenomeWideSNP_6
Total file size: 401.94 MB
RAM: 0.07MB
Number of arrays: 56
Number of references: <all arrays> (100%)
Additional parameters: [1] {flavor: chr "v2"}

dsCList <- process(cmt, verbose=verbose);
print(dsCList);

$total
AromaUnitTotalCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2
Number of files: 56
Names: A1761-01, A1761-02rh, A1761-03, ..., A1974-10 [56]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2/GenomeWideSNP_6
Total file size: 401.94 MB
RAM: 0.07MB

$fracB
AromaUnitFracBCnBinarySet:
  Name: responsify
Tags: ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2
Full name: responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2
Number of files: 56
Names: A1761-01, A1761-02rh, A1761-03, ..., A1974-10 [56]
Path (to the first file): totalAndFracBData/responsify,ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2/GenomeWideSNP_6
Total file size: 401.94 MB
RAM: 0.07MB


################################################################
################################################################
#
# Segmentation with library PSCBS
# 
# 
################################################################
################################################################
library(PSCBS)
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

#extract (total, beta) for samples
dataC = extractPSCNArray(dsCList$total)

#Segmentation all samples
pdf(paste(prefix.out, 'BetaT.pdf', sep=''), paper = 'a4') 
par(mfrow = c(2,1))
for (arr in 1:dim(dataC)[3]){
  
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  betaT = dataC[,'fracB',samp]
  plotDensity(betaT, xlim = c(0, 1), main = arr)
  grid()
}  
dev.off()

tauAs = numeric(56)
tauAs[c(1, 3, 5, 7, 9, 10, 12, 13, 18, 19:23, 27:28, 31:38, 40:42,48, 50:51, 54)] = .2
tauAs[c(2, 4, 6, 8, 13, 14, 16, 17, 24, 29:30, 39, 46:47, 49, 52:53)] = .25
tauAs[c(11, 26, 52)] = .3
tauAs[c( 15, 25, 43:45, 55:56)] = .15


yR <- getAverageFile(dsCList$total);
for (arr in 1:dim(dataC)[3]){
  
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  yT = dataC[,'total',samp]
  CT = 2 * yT/yR[,1]
  betaT = dataC[,'fracB',samp]
  ugp = getAromaUgpFile(dsNList$total)
  chromosome = ugp[,1,drop = T]
  x = ugp[, 2, drop = T]
  df = data.frame(chromosome = chromosome, x=x, CT=CT[,1], betaT = betaT)
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
  
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  saveObject(fit, file=filename);
  png(file=paste(prefix.out, 'CBSfits/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
  rm(list = c('fit','kappa','deltaAB','deltaLOH','deltaCN'))
}  


save.image(paste(prefix.out, 'SNP.RData', sep = '')) 

load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

##NOTE: NTCN calls are not optimal, so do not use them!
##Print segmentations without the NTCN cutoffs:
for (arr in 1:dim(dataC)[3]){
  
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  fit$output$ntcnCall = NULL
  png(file=paste(prefix.out, 'CBSfits/CNsegments_samp', arr, 'noNTCN.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
} 

##Number of segments in each sample
nsegs = NULL
for (arr in 1:dim(dataC)[3]){
  #Prepare data for this sample
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  nsegs = c(nsegs, nrow(fit$output))
} 
names(nsegs) = 1:56

nsegs
1     2     3     4     5     6     7     8     9    10    11    12    13 
591   819   534   702  1219   363   704   602   688   646   509   640   686 
14    15    16    17    18    19    20    21    22    23    24    25    26 
639  1067   497   695   867   544   358   875   567 32378  1016   885  1196 
27    28    29    30    31    32    33    34    35    36    37    38    39 
1126   738   496   821   545   709  1129   713   671   619   759   760 34118 
40    41    42    43    44    45    46    47    48    49    50    51    52 
519   748   554   773   619   964   548   790   644   783   709   896   858 
53    54    55    56 
581   602  1004   702 


################################################################
################################################################
#
# Data for integer copy numbers from TAPS
# 
# 
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/output/SNP/", sep='')
prefix.taps = paste(getwd(), 'TAPS/', sep = '/')


load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

### Array names
###
filenames = list.files(prefix.raw)
nc = nchar(filenames, type = 'c')
filenames = filenames[substr(filenames, nc-2, nc) == 'CEL']
filenames = substr(filenames, 1, nchar(filenames, type = 'c')-4)
info = data.frame(arrname = filenames)
info$AROS.ID2 = substr(info$arrname, 1, 8)
tmp = read.delim(paste(prefix.ann, 'A1761_work sheet_ SNP 6.0_okay_Scan_info.txt', sep = ''),
                 sep = '\t', stringsAsFactors = F)
tmp$AROS.ID2 = unlist(lapply(strsplit(tmp$AROS.ID, split = '-', fixed = T), paste, collapse = '.'))
info = merge(info, tmp, all.x = T)[,c('AROS.ID','AROS.ID2', 'arrname','Sample.ID')]
info$Sample.ID[info$AROS.ID2 %in% c('A1974.08','A1974.09','A1974.10')] = 
  c('3860 T','4541 T','4868 T')
info$AROS.ID[info$AROS.ID2 %in% c('A1974.08','A1974.09','A1974.10')] = 
  c('A1974-08','A1974-09','A1974-10')
info$arrname2 = unlist(lapply(strsplit(as.character(info$arrname), 
                                       split = '.', fixed = T), paste, collapse = '-'))
info$shortname = 1:56
tmp = info[,c('shortname','Sample.ID','arrname2')]
names(tmp) = c('Number','Sample.ID','SNP.array')
rownames(tmp) = NULL
write.table(tmp, file = paste(prefix.out, '../Key.txt', sep = ''), quote = F, sep = '\t', row.names = F) 
write.table(info, file = paste(prefix.ann, '../../Infokey.txt', sep = ''), quote = F, sep = '\t', row.names = F) 


# Identification of units in Chr 2:81-86MB and their positions
cdf <- getCdf(cesN) 
gi <- getGenomeInformation(cdf)
pos <- getAromaUgpFile(gi)


for (arr in setdiff(1:56), c(23, 39)){
     #Samples 23 and 39 were run with Nexus segments instead
  #Prepare data for this sample
  samp = as.character((info$arrname2)[arr])
  
  yT = dataC[,'total',samp]
  log2 = log2( yT/yR[,1]) #Warnings OK
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
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
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

#This was used temporarily to copy a few example files into a separate dir
prefix.tapsex = paste(getwd(), 'TAPSex/', sep = '/')
for (arr in 1:56){
  if (!(arr %in% c(3, 10, 14, 16, 22, 28, 30, 36, 54, 39, 23, 26))){
    samp = as.character((info$arrname2)[arr])
    sample.name <- paste(arr, (substr(as.character(info$Sample.ID), 1, 4))[arr], sep = '.')
    mydir = paste(prefix.taps, sample.name, sep ='')
    file.copy(paste(mydir, '/', sample.name,'_karyotype.chr1.jpeg', sep = ''),
              paste(prefix.tapsex, sample.name,'_karyotype.chr1.jpeg', sep = ''))
  }
} 
prefix.tapsex = paste(getwd(), 'TAPSex/', sep = '/')
for (arr in 1:29){
  if (!(arr %in% c(3, 4, 10, 14, 22, 28, 30, 36, 54, 6))){
    samp = as.character((info$arrname2)[arr])
    sample.name <- paste(arr, (substr(as.character(info$Sample.ID), 1, 4))[arr], sep = '.')
    mydir = paste(prefix.taps, sample.name, sep ='')
    for (chrom in c(as.character(1:22), 'X')){
      file.copy(paste(mydir, '/', sample.name,'_karyotype.chr',chrom,'.jpeg', sep = ''),
              paste(prefix.tapsex, sample.name,'_karyotype.chr', chrom, '.jpeg', sep = ''))
    }
    file.copy(paste(mydir, '/', sample.name,'_overview.jpeg', sep = ''),
              paste(prefix.tapsex, sample.name,'_overview.jpeg', sep = ''))
  }
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
prefix.taps = paste(getwd(), 'TAPS/', sep = '/')
TAPS_plot(directory = prefix.taps)

#To do TAPS_call:
#Rename SampleData.txt in last TAPS sample directory to sampleData.txt
#and put it directly under the TAPS directory!
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

myTAPS_call(directory = prefix.taps)

################################################################
################################################################
#
# For GISTIC
# 
# 
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/output/SNP/", sep='')


################################################################
# Initial check of input data format

library(PSCBS)
arr = 56  #Prepare data for this sample
fitname = paste('fit', arr, sep = '')
filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
fit <- loadObject(filename);
segtmp = getSegments(fit)

x = log2(segtmp$tcnMean)-1
smoothScatter(x, main = 'log2(2*R)-1')
grid()

x = log2(segtmp$tcnMean/2)-1
smoothScatter(x, main = 'log2(R)-1', ylim = c(-4,4))
grid()

ex = read.delim(paste(prefix.raw, '../../gisticExample/segmented_data_080520.seg', sep = ''), sep = '\t')
smoothScatter(ex$Seg.CN, main = 'ex')
grid()


################################################################
# Produce input files for GISTIC


### Markers file
###
# Identification of units and their positions
load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma)
date = format(Sys.Date())

cdf <- getCdf(cesN) 
gi <- getGenomeInformation(cdf)
pos <- getAromaUgpFile(gi)
un = getUnitNames(cdf)
out = data.frame("Marker Name" = un, Chromosome = pos[,1], 
                 "Marker Position" = pos[,2], check.names = F)
names(out)[2:3] = c('Chromosome','Marker Position')
out = out[apply(!is.na(out), 1, all),] #Delete all rows with any missing value
out = out[order(out$Chromosome, out[,'Marker Position']),]
write.table(out, paste(prefix.raw, '../../../GISTIC/Markers_file.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')
#out = read.delim(paste(prefix.raw, '../../../GISTIC/Markers_file.txt', sep = ''))
#names(out) = c('Marker Name','Chromosome','Marker Position')

### Segmentation file
###
library(PSCBS)
info$ntcn = T
info$ntcn[info$pass == F] = NA
segs = NULL
for (arr in 1:56){
  #Prepare data for this sample
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  segtmp = getSegments(fit)
  sample.name <- paste(arr, (substr(as.character(info$Sample.ID), 1, 4)[info$pass])[arr], sep = '.')
  segtmp$Sample = sample.name
  if (is.null(segtmp$ntcnCall)) {
    segtmp$ntcnCall = NA
    info$ntcn[info$Sample.ID == paste((substr(as.character(info$Sample.ID), 1, 4)[info$pass])[arr], 'T')] = F
  }
  segs = rbind(segs, segtmp)
}

out1 = data.frame(Sample = segs$Sample, Chromosome = segs$chromosome, "Start Position" = segs$tcnStart,
                  "End Position" = segs$tcnEnd, "Num markers" = segs$tcnNbrOfSNPs, 
                  Seg.CN = log2(segs$tcnMean)-1, check.names = F)

out1 = out1[apply(!is.na(out1), 1, all),] #Delete all rows with any missing value
out1 = out1[order(out1$Chromosome, out1[,'Start Position']),]
out1 = subset(out1, out1[,'Num markers'] >0)
out1 = subset(out1, !(Seg.CN %in% c(Inf, -Inf)))

lastPos = function(x, y){
  new = x
  for (i in 1:length(x)) {
    if (!x[i] %in% y) new[i] = max(y[y<x[i]])
  }
  new
}
firstPos = function(x, y){
  new = x
  for (i in 1:length(x)) {
    if (!x[i] %in% y) new[i] = min(y[y>x[i]])
  }
  new
}
for (i in 1:nrow(out1)){
  chr = out1$Chromosome[i]
  out1[i, 'Start Position'] = firstPos(out1[i, 'Start Position'], 
                                      out[out$Chromosome == chr, 'Marker Position'])
  out1[i, 'End Position'] = lastPos(out1[i, 'End Position'], 
                                       out[out$Chromosome == chr, 'Marker Position'])
}
write.table(out1, paste(prefix.raw, '../../../GISTIC/Segmentation_file.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')

#out1 = read.delim(paste(prefix.raw, '../../../GISTIC/Segmentation_file.txt', sep = ''))
save.image(paste(getwd(), '/GISTIC/GISTIC.RData', sep = ''))
load(paste(getwd(), '/GISTIC/GISTIC.RData', sep = ''))

table(out1$Chromosome)
out1$Chromosome = as.character(out1$Chromosome)
out1$Chromosome[out1$Chromosome == '23'] = 'X'
out1 = subset(out1, !(Chromosome %in% c('24','25')))
write.table(out1, paste(prefix.raw, '../../../GISTIC/Segmentation2_file.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')
out$Chromosome = as.character(out$Chromosome)
out$Chromosome[out$Chromosome == '23'] = 'X'
out = subset(out, !(Chromosome %in% c('24','25')))
write.table(out, paste(prefix.raw, '../../../GISTIC/Markers_file2.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')
################################################################
################################################################
#
# Run GISTIC
# 
# 
################################################################
################################################################

#FOR EAMPLE
#!/bin/bash
export thisdir=`pwd`

# --- RUN GISTIC 2.0 ---
echo --- creating output directory ---
  export basedir=$thisdir/example_results_IL
mkdir -p $basedir 

echo --- running GISTIC ---
export segfile=/usr/local/bioinfsoftware/gistic/current/examplefiles/segmentationfile.txt
export markersfile=/usr/local/bioinfsoftware/gistic/current/examplefiles/markersfile.txt
export refgenefile=/usr/local/bioinfsoftware/gistic/current/refgenefiles/hg16.mat
export alf=/usr/local/bioinfsoftware/gistic/current/examplefiles/arraylistfile.txt
export cnvfile=/usr/local/bioinfsoftware/gistic/current/examplefiles/cnvfile.txt

/usr/local/bioinfsoftware/gistic/current/gp_gistic2_from_seg -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -alf $alf -cnv $cnvfile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90


#FOR MY SAMPLES, run in GISTIC dir
export thisdir=`pwd`

# --- RUN GISTIC 2.0 ---
echo --- creating output directory ---
  export basedir=$thisdir/results
mkdir -p $basedir 

echo --- running GISTIC ---
  export segfile=$thisdir/Segmentation_file.txt
export markersfile=$thisdir/Markers_file.txt
export refgenefile=/usr/local/bioinfsoftware/gistic/current/refgenefiles/hg19.mat

/usr/local/bioinfsoftware/gistic/current/gp_gistic2_from_seg -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -smallmem 0 -broad 1 -conf 0.90 -rx 0 -twosides 1 -savegene 1 -savedata 0 -v 30 > log.txt

##Online instead: Proc 719480
GISTIC_2.0.result <- run.analysis(gp, "urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00125:6", refgene.file="hg19_with_miR_20120227.mat", seg.file="Segmentation2_file.txt", markers.file="Markers_file2.txt", array.list.file="", cnv.file="", gene.gistic="1", amplifications.threshold="0.1", deletions.threshold="0.1", join.segment.size="4", qv.thresh="0.25", remove.X="0", confidence.level="0.9", run.broad.analysis="0", broad.length.cutoff="0.98", max.sample.segs="2500", arm.peel="0", output.prefix="WithXp0.9")

################################################################
################################################################
#
# For ABSOLUTE with TCN
# 
# 
################################################################
################################################################
library(aroma.affymetrix)
library(PSCBS)
prefix.ABS = paste(getwd(), '/ABSOLUTE/RunAbsoluteTCN/', sep = '')
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname


for (i in 1:56){
  sample.name = paste(i, info$sample[i], sep = '.')
  fitname = paste('fit', i, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  seg.dat.fn = paste(prefix.ABS, 'seg.dat.files/', fitname, 'tcn.txt', sep = '');
  if (!file.exists(seg.dat.fn)){
    fit = loadObject(filename)
    segs = getSegments(fit)
    out = data.frame(Chromosome = segs$chromosome, Start = segs$tcnStart, End = segs$tcnEnd, 
                     Num_Probes = segs$tcnNbrOfSNPs, Segment_Mean = log2(segs$tcnMean/2))
    write.table(out, seg.dat.fn, 
                row.names = F, quote = F, sep = '\t')
    
  }
}

################################################################
################################################################
#
# For ASCAT
# 
# 
################################################################
################################################################

setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/output/SNP/", sep='')
prefix.taps = paste(getwd(), 'TAPS/', sep = '/')


load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))

# Identification of units and their positions
cdf <- getCdf(cesN) 
gi <- getGenomeInformation(cdf)
pos <- getAromaUgpFile(gi)
unf <- getUnitNamesFile(pos)
unitNames <- getUnitNames(unf)

out = data.frame(chrs = pos[,1], pos = pos[,2])
names(out) = c('chrs','pos')
out$chrs = as.character(out$chrs)
out$chrs[out$chrs == '23'] = 'X'
out$chrs[out$chrs == '24'] = 'Y'
out$chrs[out$chrs == '25'] = 'M'
rownames(out) = unitNames   
for (arr in 1:56){
  #Prepare data for this sample
  samp = as.character((info$arrname2)[arr])
  yT = dataC[,'total',samp]
  log2 = log2( yT/yR[,1]) #Warnings OK
  cname = paste('S', arr, sep = '')
  out = cbind(out, log2)
  colnames(out)[ncol(out)] = cname
}  
d = file.path(getwd(), 'ASCAT', 'Indata')
write.table(out, file.path(d, 'Tumour_lrr.txt'), quote = F,
            sep = '\t')


out = data.frame(chrs = pos[,1], pos = pos[,2])
names(out) = c('chrs','pos')
out$chrs = as.character(out$chrs)
out$chrs[out$chrs == '23'] = 'X'
out$chrs[out$chrs == '24'] = 'Y'
out$chrs[out$chrs == '25'] = 'M'
rownames(out) = unitNames   
for (arr in 1:56){
  #Prepare data for this sample
  samp = as.character((info$arrname2)[arr])
  betaT = dataC[,'fracB',samp]
  cname = paste('S', arr, sep = '')
  out = cbind(out, betaT)
  colnames(out)[ncol(out)] = cname
}  
d = file.path(getwd(), 'ASCAT', 'Indata')
write.table(out, file.path(d, 'Tumour_baf.txt'), quote = F,
            sep = '\t')






################################################################
################################################################
#
# ASCAT segmentation
# 
# 
################################################################
################################################################

d = file.path(getwd(), 'ASCAT', 'ASCAT2.1')
source(file.path(d, "ascat.R"))
d = file.path(getwd(), 'ASCAT', 'Indata')
ascat.bc <- ascat.loadData(Tumor_LogR_file=file.path(d, 'Tumour_lrr.txt'), 
                           Tumor_BAF_file=file.path(d, 'Tumour_baf.txt'))

d = file.path(getwd(), 'ASCAT', 'ASCAT2.1')
source(file.path(d,"predictGG.R"))
save.image('ASCAT/Ascat.RData')
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")
save.image('Ascat_gg.RData')
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)
save.image('Ascat_bc.RData')
ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc)
save.image('Ascat_output.RData')


source(file.path(d,"organize.ascat.output-1.R"))
tumor.segments <- organize.ascat.segments(ascat.output, ascat.bc$SNPpos)
write.table(tumor.segments,file="./Tumor_segments.txt",quote=FALSE,sep="\t")
print("Tumor analysis completed")
save(ascat.output = ascat.output, file = 'ascat.output.RData')

###


######################################################################################
######################################################################################
