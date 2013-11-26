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
library(crlmm)
library(aroma.affymetrix)

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
  fit <- segmentByNonPairedPSCBS(df, knownSegments=knownSegments,
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
  } else fit = callB(fit, verbose = -100)
  
  ## Call loss of heterozygosity (LOH)
  ## (a) Estimate C_1 threshold for calling LOH
  deltaLOH <- estimateDeltaLOH(fit); # If skipped, will be done internally
  
  ## (b) Call LOH based on bootstrapped segment C_1 levels
  fit <- callLOH(fit, delta=deltaLOH, verbose=-10);
  
  ## Call NTCN
  ## (a) Estimate the threshold for calling neutral TCN segments
  ##     By shrinking 'scale', more segments will be non-NTCN.
  deltaCN <- estimateDeltaCN(fit, scale=1.0, kappa=kappa);
  
  ## (b) Call NTCN based on bootstrapped segment TCN levels
  fit <- callNTCN(fit, delta=deltaCN, verbose=-20);
  
  #  toPNG(fitname, tags="tracks,avgDH=median,AB+LOH+NTCN", aspectRatio=0.6, {
  #    plotTracks(fit, chromosomes=1:23);
  #  });
  
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  saveObject(fit, file=filename);
  png(file=paste(prefix.out, 'CBSfits/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
  
  
  
  
  
  
  fit = segmentByNonPairedPSCBS(df, knownSegments = knownSegments, avgDH="median", 
                                                     tauA = tauAs[arr], seed = 48879, verbose = -10) 

  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  saveObject(fit, file=filename);
  png(file=paste(prefix.out, 'CBSfits/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
}
save.image(paste(prefix.out, 'SNP.RData', sep = '')) 

#Call NTCN segments
for (arr in 1:dim(dataN)[3]){
  #Prepare data for this sample
  samp = as.character((info$shortname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  
  deltaCN <- estimateDeltaCN(fit, scale=1) #0.2601609 for arr = 1
  fit1 = try(callAB(fit, delta = deltaCN, verbose = -10), silent = T)
  if (substr(fit1[1], 1, 5) != 'Error') fit = fit1
  fit1 = try(callNTCN(fit, verbose = -10), silent = T)
  if (substr(fit1[1], 1, 5) != 'Error') fit = fit1
  saveObject(fit, file=filename);
  seg = getSegments(fit, simplify = T)
  saveObject(seg=seg, file = paste(prefix.out, 'CBSfits/CNsegments_samp', arr, 
                               '.RData', sep = ''))
}

save.image(paste(prefix.out, 'SNP.RData', sep = '')) 

load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())

library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

sampleName <- "ingrid";

## Load the data from a previously saved segmentation object.
pathname <- sprintf("%s.Rbin", sampleName);
fit <- loadObject(pathname);
data <- getLocusData(fit);
rm(fit);


## Locate centromeres and other large gaps in data
gaps <- findLargeGaps(data, minLength=1e6);
knownSegments <- gapsToSegments(gaps);
## => dim(knownSegments) == c(64, 4)

## 'betaT' is already normalized using TumorBoost => tbn=FALSE.
fit <- segmentByNonPairedPSCBS(data, knownSegments=knownSegments,
                               avgDH="median", seed=0xBEEF, verbose=-10);
printf("Number of segments: %d\n", nbrOfSegments(fit));
## Number of segments: 568

toPNG(sampleName, tags="tracks,avgDH=median", aspectRatio=0.6, {
  plotTracks(fit);
});


## Estimate global background level in [0,1] (due to normal
contamination and more)
kappa <- estimateKappa(fit, verbose=-10);
printf("Kappa: %g\n", kappa);
## Kappa: 0.374814


## Call allelic balance
## (a) Estimate DH threshold for calling AB
deltaAB <- estimateDeltaAB(fit, kappa=kappa); # If skipped, will be done internally
printf("Delta_AB: %g\n", deltaAB);
## Delta_AB: 0.146787

## (b) Call AB based on bootstrapped segment DH levels
fit <- callAB(fit, delta=deltaAB, verbose=-100);
printf("Number of AB segments: %d\n", sum(getSegments(fit)$abCall, na.rm=TRUE));
## Number of AB segments: 193

toPNG(sampleName, tags="tracks,avgDH=median,AB", aspectRatio=0.6, {
  plotTracks(fit);
});


## Call loss of heterozygosity (LOH)
## (a) Estimate C_1 threshold for calling LOH
deltaLOH <- estimateDeltaLOH(fit); # If skipped, will be done internally
printf("Delta_LOH: %g\n", deltaLOH);
## Delta_LOH: 0.566545

## (b) Call LOH based on bootstrapped segment C_1 levels
fit <- callLOH(fit, delta=deltaLOH, verbose=-10);
printf("Number of LOH segments: %d\n", sum(getSegments(fit)$lohCall,
                                           na.rm=TRUE));
## Number of LOH segments: 149

toPNG(sampleName, tags="tracks,avgDH=median,AB+LOH", aspectRatio=0.6, {
  plotTracks(fit);
});


## Call NTCN
## (a) Estimate the threshold for calling neutral TCN segments
##     By shrinking 'scale', more segments will be non-NTCN.
deltaCN <- estimateDeltaCN(fit, scale=1.0, kappa=kappa);
printf("Delta_CN: %g\n", deltaCN);
## Delta_CN: 0.312593

## (b) Call NTCN based on bootstrapped segment TCN levels
fit <- callNTCN(fit, delta=deltaCN, verbose=-20);
printf("Number of NTCN segments: %d\n", sum(getSegments(fit)$ntcnCall,
                                            na.rm=TRUE));
## Number of NTCN segments: 303

toPNG(sampleName, tags="tracks,avgDH=median,AB+LOH+NTCN", aspectRatio=0.6, {
  plotTracks(fit);
});


#### Printed this for preliminary overview
####
load(paste(prefix.out, 'SNP.RData', sep = ''))
library(PSCBS)
date = format(Sys.Date())

for (arr in 1:dim(dataN)[3]){
    #Prepare data for this sample
    fit <- get(paste('fit', arr, sep = ''))
    seg = getSegments(fit, simplify = T)
    save(seg=seg, file = paste(prefix.out, 'CBSfits/NoCalmate/CNsegments_samp', arr, 
                           '.RData', sep = ''))
}







################################################################
################################################################
#
# Segmentation with library PSCBS: new attempt
# 
# 
################################################################
################################################################
library(PSCBS)
library("R.devices");
devOptions("png", width=1024);
setOption("devEval/args/force", FALSE);

for (arr in 5:56){
  
  #Prepare data for this sample
  samp = as.character((info$shortname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  df <- getLocusData(fit);
  rm(fit);
  
  gaps = findLargeGaps(df, minLength = 1e+06)
  knownSegments = gapsToSegments(gaps)
  
  ## 'betaT' is already normalized using TumorBoost => tbn=FALSE.
  fit <- segmentByNonPairedPSCBS(df, knownSegments=knownSegments,
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
  } else fit = callB(fit, verbose = -100)
  
  ## Call loss of heterozygosity (LOH)
  ## (a) Estimate C_1 threshold for calling LOH
  deltaLOH <- estimateDeltaLOH(fit); # If skipped, will be done internally
  
  ## (b) Call LOH based on bootstrapped segment C_1 levels
  fit <- callLOH(fit, delta=deltaLOH, verbose=-10);
  
  ## Call NTCN
  ## (a) Estimate the threshold for calling neutral TCN segments
  ##     By shrinking 'scale', more segments will be non-NTCN.
  deltaCN <- estimateDeltaCN(fit, scale=1.0, kappa=kappa);
  
  ## (b) Call NTCN based on bootstrapped segment TCN levels
  fit <- callNTCN(fit, delta=deltaCN, verbose=-20);
  
#  toPNG(fitname, tags="tracks,avgDH=median,AB+LOH+NTCN", aspectRatio=0.6, {
#    plotTracks(fit, chromosomes=1:23);
#  });
  
  filename <- paste(prefix.out, 'CBSfits/withcalls/', fitname, '.RData', sep = '');
  saveObject(fit, file=filename);
  png(file=paste(prefix.out, 'CBSfits/withcalls/CNsegments_samp', arr, '.png', sep = ''), width=800, height=580)
  plotTracks(fit, chromosomes=1:23)
  title(main = samp, outer = T)
  dev.off()   
}

save.image(paste(prefix.out, 'SNPwithcalls.RData', sep = '')) 

load(paste(prefix.out, 'SNP.RData', sep = ''))
library(aroma.affymetrix)
library(PSCBS)
date = format(Sys.Date())


################################################################
################################################################
#
# CRLMM
#
################################################################
################################################################


library(ff)
library(genefilter)
library(IRanges)
library(MASS)
library(VanillaICE)
library(crlmmCompendium)
outdir = prefix.dat
ldPath(prefix.out)
ocProbesets(50000)
ocSamples(200)
library(cacheSweave)
setCacheDir(outdir)


######### CRLMM
#########

crlmmResults = crlmm(info$filenames[info$pass == T], gender = info$Gender[info$pass == T])


######### Summarize results with more info
#########

dim(calls(crlmmResults))
[1] 906600     56

calls(crlmmResults)[1:10, 1:3]
confs(crlmmResults, FALSE)[1:10, 1:3]

i2p(confs(crlmmResults, FALSE)[1:10, 1:3]) #P-values are huge, so genotypes are not very useful.





######################################################################################
######################################################################################
#END OF FILE 










################################################################
################################################################
#
# Aroma copy numbers step 5: (JUST JUNK EXAMPLES)
# Calculation of raw copy numbers
# This section is just examples in case I need later
################################################################
################################################################

#The above 'cesN' object contains chip-effect estimates according to the 
#CRMA v2 method (Bengtsson et al. 2008c).  
#In this section we will show how to calculate raw copy numbers relative 
#to a reference.  Note that several of the downstream methods, such as segmentation 
#methods, will do this automatically/internally and therefore it is often not necessary 
#to do this manually.

#Deciding on a reference

#There are two common use cases for CN analysis; either one do (i) a paired analysis where 
#each sample is coupled with a unique reference (e.g. tumor/normal) or (ii) a non-paired analysis 
#where each sample use the same common reference.  When a common reference is used, it is often 
#the average of a pool of samples.  Here we will show how to do the latter.

#Calculating the robust average of all samples

#To calculate the robust average of chip effects across all existing samples (i=1,2,...,I), 
#that is, theta_Rj = median_i {theta_ij} for each unit j=1,2,...,J.  
#This calculation can be done as:

ceR <- getAverageFile(cesN, verbose=verbose)
print(ceR)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Pasted onto here in unix301 run

#From the output we learn that 'ceR' is a CnChipEffectFile just like the other arrays in 
#the 'cesN' set.  From now on we can treat it as if it was the output a hybridization although 
#it is actually the average over many.  The main difference is that this one is likely to have 
#more precise chip effects (because of the averaging over many estimates).

#Extracting chip-effect estimates

#The next step is to calculate the relative copy numbers:
#C_ij = 2*theta_ij / theta_Rj
#where we assume the copy-neutral state has two (2) copies.  In order to calculate this, 
#we need to extract theta_ij and theta_Rj.

#Example: A single unit in one sample
#Consider array i=3 and unit j=987. We can then extract these values as:

ce <- getFile(cesN, 3)  # Array #3
theta <- extractTheta(ce, unit=987)
thetaR <- extractTheta(ceR, unit=987)
C <- 2*theta/thetaR
print(C)
[,1]
[1,] 1.917

#Thus, for array #3 (NA06993) and unit #987 (SNP_A-4268099), the estimated total raw copy number 
#is C=1.92.  If the DNA extract is homogeneous (containing only the same normal cells), we expect 
#this to correspond to a truly diploid locus.

#If we look at the individual theta and thetaR estimates,



print(c(theta, thetaR))
[1] 2570 2681

#we find that, on average (robust), an array has chip-effect signal 2681 (on the intensity scale).  
#This is the signal we expect a copy neutral event at this particular locus to have.

#Example: A small region on chromosome 2 in one sample
#Next we are interested in the distribution of raw copy number {C_ij} in a small region on 
#chromsome 2 in sample NA06991.

# Identification of units in Chr 2:81-86MB and their positions
cdf <- getCdf(cesN) 
gi <- getGenomeInformation(cdf)
units <- getUnitsOnChromosome(gi, chromosome=2, region=c(81,86)*1e6)
pos <- getPositions(gi, units=units)

# Retrieving CN estimates of the reference in the above region
ceR <- getAverageFile(cesN)
thetaR <- extractTheta(ceR, units=units)

# Retrieving the corresponding estimates of sample A1761-56
ce <- getFile(cesN, indexOf(cesN, "A1761-56"))
theta <- extractTheta(ce, units=units)

# Calculate the relative CNs
C <- 2*theta/thetaR

# Plotting data along genome 
filename <- sprintf("%s,Chr02,81-86Mb.png", getName(ce))
png(file=paste(prefix.out,filename, sep = ''), width=640, height=640)
par(mar=c(3,4,2,1)+0.1, mfrow=c(2,1))
plot(pos/1e6, C[,1], pch=".", cex=3, ylim=c(0,4))
stext(side=3, pos=0, getName(ce))
stext(side=3, pos=1, "Chr2")
stext(side=3, pos=.5, "minor")
plot(pos/1e6, C[,2], pch=".", cex=3, ylim=c(0,4))
stext(side=3, pos=0, getName(ce))
stext(side=3, pos=1, "Chr2")
stext(side=3, pos=.5, "major")
dev.off()
#dev.print(png, file=filename, width=640, height=320)

#Figure: Copy number estimates in a 5.0Mb region on Chr2 for sample A1761-56.  
#There is a clear deletion at 83.1-83.7Mb. 


######################################################################################
######################################################################################
