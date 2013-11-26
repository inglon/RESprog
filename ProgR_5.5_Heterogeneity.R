################################################################
################################################################
# Script: ProgR_5.5_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  23/05/2013
# R version: R 3.0.0
# Details: Derive basic heterogeneity estimates
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

prefix.raw = paste(getwd(), "/AROMA/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/AROMA/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/ABSOLUTE/HAPSEG/", sep='')


date = format(Sys.Date())

source(file.path(getwd(), 'programs', 'ProgR_0.1_Functions.R'))



################################################################
################################################################
#
# Format intensities for HAPSEG
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample = substr(info$Sample.ID, 1, 4)
info$sample.name = paste(rownames(info), info$sample, sep = '.')
info$myarray = info$arrname

setwd(paste(getwd(), '/AROMA', sep=''))

     
######Find allele intensities

for (i in 1:56){   
library(aroma.affymetrix)
dataSet <- "responsify";
tags <- "ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2";
chipType <- "GenomeWideSNP_6";
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
df <- getFile(ds, which(getNames(ds) == info$arrname[i]));
bs <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
bf <- getFile(bs, which(getNames(ds) == info$arrname[i]));
ugp <- getAromaUgpFile(df)
unf <- getUnitNamesFile(ugp)
unitNames <- getUnitNames(unf)
chromosome <- ugp[, 1, drop = TRUE]
position <- ugp[, 2, drop = TRUE]
data = data.frame(unitNames = unitNames, chromosome=chromosome, position = position, 
                  total = df[,1], fracB = bf[,1], stringsAsFactors = F)
names(data)[4:5] = c('total','fracB')
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
dir.create(d, recursive= T, showWarnings = F)
write.table(data, file.path(d, 'Total,fracB.txt'), row.names = F, quote = F, sep = '\t')
}


######Output channel specific intensities

#First move from old to new directory, only run this once!
for (i in 1:56){
  d.old = file.path(getwd(), 'ABSOLUTE', 'Old', 'HAPSEG', info$sample.name[i])
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  dir.create(d, recursive = T, showWarnings = F)
  file.copy(file.path(d.old, 'Total,fracB.txt'), file.path(d, 'Total,fracB.txt'))
  file.copy(file.path(d.old, 'Hetsegs.txt'), file.path(d, 'Hetsegs.txt'))
  file.copy(file.path(d.old, 'Homsegs.txt'), file.path(d, 'Homsegs.txt'))
  file.copy(file.path(d.old, 'Subsegs.txt'), file.path(d, 'Subsegs.txt'))
  
}

#Now output the new intensities
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  data = read.delim(file.path(d, 'Total,fracB.txt'), stringsAsFactors = F)
  
  ###Write to file for HAPSEG
  data$baf = pmax(data$fracB, 0)
  data$baf = pmin(data$fracB, 1)
  data$A = data$total*(1-data$baf)
  data$B = data$total*(data$baf)
  
  dat = subset(data, select = c('A','B'))
  rownames(dat) = data$unitNames
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  dir.create(d, recursive = T, showWarnings = F)
  save(dat, file = file.path(d, 'snpfile.RData'))
}

#Intensities for nexus preprocessed samples
#Need the previous section to just have been run, becuase will build on the object 'data'
for (i in c(23, 39)){
  d = file.path(getwd(), 'TAPS', info$sample.name[i])
  probes = read.delim(file.path(d, 'probes.txt'))
  bafs = read.delim(file.path(d, 'snps.txt'))
  data$total = NA
  data$fracB = NA
  probes$Chromosome = as.character(probes$Chromosome)
  nc = nchar(probes$Chromosome, type = 'c')
  probes$Chromosome = substr(probes$Chromosome, 4, nc)
  probes$Chromosome[probes$Chromosome == 'X'] = '23'
  probes$Chromosome[probes$Chromosome == 'Y'] = '24'
  probes$Chromosome = factor(probes$Chromosome, levels = as.character(1:25))
  bafs$Chromosome = as.character(bafs$Chromosome)
  nc = nchar(bafs$Chromosome, type = 'c')
  bafs$Chromosome = substr(bafs$Chromosome, 4, nc)
  bafs$Chromosome[bafs$Chromosome == 'X'] = '23'
  bafs$Chromosome[bafs$Chromosome == 'Y'] = '24'
  bafs$Chromosome = factor(bafs$Chromosome, levels = as.character(1:25))
  data$chromosome = factor(data$chromosome, levels = as.character(1:25))
  
  h = RangedData(IRanges(start = probes[,'Start'], width = 1),  
                 space = probes$Chromosome, universe = 'hg19', val = probes[,'Value']) 
  ix = which(!is.na(data$chromosome))
  d = RangedData(IRanges(start = data$position[ix], width = 1),  
                 space = data$chromosome[ix], universe = 'hg19')
  hits = findOverlaps(query = h, subject = d)
  
  for (j in 1:length(levels(data$chromosome))){
    ch = levels(data$chromosome)[j]
    ind = which(data$chromosome == ch & !is.na(data$chromosome))
    data$total[ind[subjectHits(hits[[ch]])]] =   (probes$Value[probes$Chromosome 
                                                                  == ch])[queryHits(hits[[ch]])]
  }
  h = RangedData(IRanges(start = bafs[,'Start'], width = 1),  
                 space = bafs$Chromosome, universe = 'hg19', val = bafs[,'Value']) 
  hits = findOverlaps(query = h, subject = d)
  
  for (j in 1:length(levels(data$chromosome))){
    ch = levels(data$chromosome)[j]
    ind = which(data$chromosome == ch & !is.na(data$chromosome))
    data$fracB[ind[subjectHits(hits[[ch]])]] =   (bafs$Value[bafs$Chromosome 
                                                               == ch])[queryHits(hits[[ch]])]
  }
  
  data$total = 2**data$total
  data$baf = pmax(data$fracB, 0)
  data$baf = pmin(data$fracB, 1)
  data$A = data$total*(1-data$baf)
  data$B = data$total*(data$baf)
  
  dat = subset(data, select = c('A','B'))
  rownames(dat) = data$unitNames
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  dir.create(d, recursive = T, showWarnings = F)
  save(dat, file = file.path(d, 'snpfile.RData'))
}



###### Segment datasets from PSCBS

library(PSCBS)
for (i in 2:56){
  print(i)
  #Prepare data for this sample
  fitname = paste('fit', i, sep = '')
  filename <- paste(getwd(), '/AROMA/reports/SNP_segmentation/CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  seg = getSegments(fit)
  seg$Sample = info$myarray[i]
  seg = data.frame(Sample = seg$Sample, Chromosome = seg$chromosome, "Start Position" = seg$tcnStart,
                   "End Position" = seg$tcnEnd, "Num markers" = seg$tcnNbrOfSNPs, 
                   check.names = F)
  seg = subset(seg, !is.na(Chromosome) & Chromosome != 25)
  seg$Chromosome = factor(seg$Chromosome, levels = as.character(1:24))
  seg$cn = NA
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i], 'Segdat.txt')
  write.table(seg, d, row.names = F, quote = F, sep = '\t')
}

#Fix the nexus normalized samples
for (i in c(23, 39)){
  d = file.path(getwd(), 'TAPS', info$sample.name[i])
  seg = read.delim(file.path(d, 'segments.txt'))
  seg = data.frame(Sample  = info$myarray[i], Chromosome = seg$Chromosome, "Start Position" = seg$Start,
                   "End Position" = seg$End, "Num markers" = NA, 
                   check.names = F)
  seg = subset(seg, !is.na(Chromosome) & Chromosome != 'chr25')
  seg$Chromosome = as.character(seg$Chromosome)
  nc = nchar(seg$Chromosome, type = 'c')
  seg$Chromosome = substr(seg$Chromosome, 4, nc)
  seg$Chromosome[seg$Chromosome == 'X'] = '23'
  seg$Chromosome[seg$Chromosome == 'Y'] = '24'
  seg$Chromosome = factor(seg$Chromosome, levels = as.character(1:24))
  seg$cn = NA
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i], 'Segdat.txt')
  write.table(seg, d, row.names = F, quote = F, sep = '\t')
  
}
################################################################
################################################################
#
# HAPSEG
#
################################################################
################################################################
library(HAPSEG)
library(foreach)
library(doMC)
registerDoMC(20)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(info$ind, info$sample, sep = '.')


DoHapseg1 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG1', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG1', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
            use.pop='CEPH', 
            impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
            normal=FALSE, 
            out.p=0.001, seg.merge.thresh=1, use.normal = FALSE, adj.atten = FALSE, 
            phased.bgl.dir, calibrate.data=FALSE, 
            drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=1:56, .combine=c) %dopar% {
  DoHapseg1(myarray = info$myarray[i], mysample = info$sample.name[i])
}

DoHapseg2 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-10, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=1:56, .combine=c) %dopar% {
  DoHapseg2(myarray = info$myarray[i], mysample = info$sample.name[i])
}

DoHapseg3 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-80, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=1:56, .combine=c) %dopar% {
  DoHapseg3(myarray = info$myarray[i], mysample = info$sample.name[i])
}




################################################################
################################################################
#
# Format exome mutation file for ABSOLUTE
#
################################################################
################################################################
# Create sample_list for Pearl script TSV2MAF_pearl_script.rtf
# and reformat mutations file to MAF (Broad's format)
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')
prefix.in = paste(getwd(), "/rawData/Exome/", sep='')
prefix.raw = paste(getwd(), "/datasets/", sep='')

###Sample_list
mut = read.delim(paste(prefix.in, 'Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                       sep = ''))
samplesT = unique(mut$SAMPLE)
samples = sort(unique(substr(as.character(samplesT), 1, 4)))
sampT = paste(samples, 'T', sep = '')
sampN = paste(samples, 'N', sep = '')
out = rbind(data.frame(sample = sampT, type = 'Tumour', pairsamp = sampN, bam = 'dummy'), 
            data.frame(sample = sampN, type = 'Normal', pairsamp = sampT, bam = 'dummy'))
write.table(out, file = paste(prefix.in, 'sample_list.txt', sep = ''), row.names = F,col.names = F,
            quote = F, sep = '\t')

###Pearl calls
###Run from Exome library
###This never worked with our system, but Jason Li run it for us.
perl parse_tsv_to_maf_format.pl sample_list.txt HUGO_EntrezIDs.txt Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt ALL_130516_Filter_CBS_MinDepth_20.maf


### Find additional variables to merge with the maf file
### AND split maf file into one file for each sample
mut = read.delim(paste(prefix.in, 'Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                       sep = ''))
dat = data.frame(Tumor_Sample_Barcode = mut$SAMPLE,
                 Chromosome = mut$X.CHROM, Start_position = mut$POS, 
                 t_alt_count = mut$X4688T_Variant_reads,
                 t_ref_count = mut$X4688T_Read_depth - mut$X4688T_Variant_reads)
dat = unique(dat)
mafbase = 'ALL_130516_Filter_CBS_MinDepth_20'
maf = read.delim(paste(prefix.in, mafbase,'.maf', sep = ''), 
                 row.names = NULL, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", "---"), 
                 blank.lines.skip = TRUE, comment.char = "#")
names(maf)[names(maf) == 'Start_Position'] = 'Start_position'
names(maf)[names(maf) == 'End_Position'] = 'End_position'
nrow(maf)
maf = merge(maf, dat, all = T, sort = F)
nrow(maf)

info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(as.character(info$Sample.ID), 1, 4)
for (i in 1:nrow(info)){
  samp = info$sample[i]
  out = subset(maf, Tumor_Sample_Barcode == paste(samp, 'T', sep = ''))
  if (nrow(out)>0) write.table(out, file = paste(prefix.in,'MAF/', 
              mafbase, '_', i,'.', samp, '.maf', sep = ''), row.names = F, quote = F,
                               sep = '\t')
}


#maf = read.delim(maf.fn, 
#                 row.names = NULL, stringsAsFactors = FALSE, 
#                 check.names = FALSE, na.strings = c("NA", "---"), 
#                 blank.lines.skip = TRUE, comment.char = "#")





################################################################
################################################################
#
# Run ABSOLUTE without scaling, just for comparison
#
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')

DoAbsolute1 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG1/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute1', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute1/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 0.95, max.ploidy = 12, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute1(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}


DoAbsolute2 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
     '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 0.95, max.ploidy = 12, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute2(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoAbsolute3 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG3/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute3', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute3/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 0.95, max.ploidy = 12, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute3(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}


#Not used yet
d = paste(getwd(), '/ABSOLUTE/RunAbsolute1', sep = '')
obj.name = 'ALL'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                    info$sample.name[c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56)], 
                                                           '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)


################################################################
################################################################
#
# Pick patterns
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')

dir.create('HAPSEG1')
for (i in 1:56){
  dir.create(file.path('HAPSEG1', info$sample.name[i]))
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG1', info$sample.name[i])
  f = paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')
  file.copy(file.path(d, f), file.path('HAPSEG1', info$sample.name[i], f))
}
dir.create('HAPSEG3')
for (i in 1:56){
  dir.create(file.path('HAPSEG3', info$sample.name[i]))
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', info$sample.name[i])
  f = paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')
  file.copy(file.path(d, f), file.path('HAPSEG3', info$sample.name[i], f))
}


for (i in 1:56){
  
  i = 15
  
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', info$sample.name[i])
  da = load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
  da = as.data.frame(seg.dat$allele.segs)
  da$Chromosome = factor(da$Chromosome, levels = as.character(1:22))
  names(da)[6:7] = c('A1','A2')
  da$tot = (da$A1+da$A2)
  da$imba = 2*abs(da$A1/(da$A1+da$A2) - 0.5)
  d0.imba = 1/200
  d0.tot = diff(range(quantile(da$tot, c(.05, .95))))/500
  n = 100
  
  out = NULL
  for (pattern in 1:20){
    plot(da$tot, da$imba, pch = 16, col = col.transp('grey', .3), cex = 1,
         xlab = 'Total Copy Number', ylab = 'Alleleic imbalance 2*|BAF - 0.5|',
         main = info$sample.name[i])#, xlim = c(0,5))
    grid()
    #Start picking a new pattern
    print('Choose points within this pattern, press Esc when finished!')
    pos = locator(n)   #Choose points within this specific pattern, press Esc when finished!
    di = 0
    dt = 0
    tmp = NULL
    plot(da$tot, da$imba, pch = 16, col = col.transp('grey', .3), cex = 1,
         xlab = 'Total Copy Number', ylab = 'Alleleic imbalance 2*|BAF - 0.5|',
         main = info$sample.name[i])
    grid()
    
    #########RUN THESE LINES UNTIL JUST TOO MANY SEGMENTS HIGHLIGHTED
    for (findpoints in 1:1000){
      if (di>0) last.tmp = tmp
      di = di + d0.imba
      dt = dt + d0.tot
      y.close = outer(da$imba, pos$y, function(x, y, di){abs(x-y)<di},
                      di = di)
      x.close = outer(da$tot, pos$x, function(x, y, dt){abs(x-y)<dt},
                      dt = dt)
      ix = apply(y.close & x.close, 1, any)
      tmp = da[ ix,]
      points(tmp$tot, tmp$imba, col = 'blue', pch = 16, cex = 1)
      points(pos$x, pos$y, pch = 4, col = 'red')
      print('Continue [Return] or give pattern name [e.g. maj2, min1, tot4 or bal]:')
      finished.this.pattern = readline()
      if (finished.this.pattern != '') {
        last.tmp$pattern = finished.this.pattern
        break
      }
    }
    
    #########THEN SAVE THE PREVIOUS SEGMENTS (RUN THIS ONCE)
    rownames(last.tmp) = NULL
    out = rbind(out, last.tmp)    
    
    print('Find more patterns? (y/n)')
    finished.all.patterns = readline()
    if (finished.all.patterns == 'n') break
  }
  
  plot(da$tot, da$imba, pch = 16, col = col.transp('grey', .3), cex = 1,
       xlab = 'Total Copy Number', ylab = 'Alleleic imbalance 2*|BAF - 0.5|',
       main = info$sample.name[i])
  grid()
  for (p in seq_along(unique(out$pattern)))  points(imba~tot, data = out[out$pattern == unique(out$pattern)[p],], 
                                                    col = col.transp(p, .5), pch = 16)
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  dev.print(png, file=file.path(d, 'Patterns.png'), 
            width=640, height=640)
  write.table(out, file.path(d, 'Patternsegs.txt'), 
              row.names = F, quote = F, sep = '\t')
}

################################################################
################################################################
#
# Pick TAPS clusters
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')



for (i in 1:56){
  
  i = 56
  ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                            paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  
  d0.imba = diff(range(ta$imba, na.rm = T))/200
  d0.log = diff(range(ta$log2, na.rm = T))/1000
  out = NULL
  
  n = 1 ###MANUALLY INSERT THE NUMBER OF CLUSTERS YOU LIKE HERE:
  pos = locator(n)
  
  ###FOR EACH CLUSTER (CHANGE c BETWEEN 1, 2, 3, ..., n)
  c = 1
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  di = 0
  dl = 0
  
  #########RUN THESE LINES UNTIL JUST TOO MANY SEGMENTS HIGHLIGHTED
  if (di>0) last.tmp = subset(ta, abs(ta$imba - pos$y[c])<di & abs(ta$log2 - pos$x[c])<dl)
  di = di + d0.imba
  dl = dl + d0.log
  tmp = subset(ta, abs(ta$imba - pos$y[c])<di & abs(ta$log2 - pos$x[c])<dl)
  points(imba~log2, data = tmp, col = c+1, pch = 16)
  
  #########THEN SAVE THE PREVIOUS SEGMENTS (RUN THIS ONCE)
  out = rbind(out, cbind(last.tmp, cluster = c))    
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i], bg = 'white')
  points(imba~log2, data = last.tmp, col = c+1, pch = 16)
  
  
  
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  points(imba~log2, data = out, col = out$cluster+1, pch = 16)
  d = file.path(getwd(), 'ABSOLUTE', 'Choose_solution')
  dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Chosen_subclones.png', sep = '')), 
            width=640, height=640)
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  dir.create(d, recursive = T, showWarnings = F)
  write.table(out[,c('Chromosome','Start','End','cluster')], file.path(d, 'Subsegs.txt'), 
              row.names = F, quote = F, sep = '\t')
}


################################################################
################################################################
#
# Plot PSCBS: Allele1 vs Allele2 coloured by the chosen clusters
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')

library(PSCBS)

for (i in 1:56){
d = file.path(getwd(), 'AROMA', 'reports', 'SNP_segmentation', 'CBSfits')
da = loadObject(paste(d, '/fit', i, '.RData', sep = ''))
seg = da$output
#plot(dhMean~log(tcnMean), data = seg,  main = info$sample.name[i], xlim = c(-.2, 2))
#plot(c1Mean~c2Mean, data = seg, main = info$sample.name[i])
grid()
mat = seg
names(mat)[names(mat) == 'c1Mean'] = 'A1'
names(mat)[names(mat) == 'c2Mean'] = 'A2'
names(mat)[names(mat) == 'tcnStart'] = 'Start.bp'
names(mat)[names(mat) == 'tcnEnd'] = 'End.bp'
names(mat)[names(mat) == 'chromosome'] = 'Chromosome'

mat$length = mat$End.bp - mat$Start.bp
mat = subset(mat, !is.na(length) & Chromosome<23)

###Read selected segments
d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
tmp = getKnownClustersAndOverlaps(mat = mat, d = d)
cl = tmp$cl
clust = tmp$clust
w = tmp$w  


#A1 vs A2
plot1v2(a1 = mat$A1, a2=mat$A2, mat = mat, clust = clust, cl = cl, w = w)
title(main = i)
plot1v2(a1 = mat$dhMean, a2=log(mat$tcnMean), mat = mat, clust = clust, cl = cl, w = w)
title(main = i)
}
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '.............png', sep = '')), 
          width=640, height=640)
#This plot, too has the problem with CN 2 clusters not on a vertical line.



################################################################
################################################################
#
# Plot HAPSEG data with chosen clusters coloured
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')

#Either plot in one large pdf or plot each sample in a separate pdf.
#d = file.path(getwd(), 'ABSOLUTE', NA)
#pdf(paste(d, 'Cluster_overview2', '.pdf', sep=''), 
#    paper = 'a4', height = 12, width = 8)  

for (i in 1:56){
  
  #Path to print output
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData')
  pdf(file.path(d, info$sample.name[i], 'HAPSEG3_choosen_clusters.pdf'), 
      paper = 'a4', height = 12, width = 8)  
  
  
  ### Load data for this sample
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', info$sample.name[i])
  f = load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
  ab = as.data.frame(seg.dat$allele.segs)
  ab$Chromosome = factor(ab$Chromosome, levels = as.character(1:22))
  names(ab)[6:7] = c('A1','A2')
  
  ###Read selected segments
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
  tmp = getKnownClustersAndOverlaps(mat = ab, d = d)
  cl = tmp$cl
  clust = tmp$clust
  w = tmp$w  
  
  ### Plots of HAPSEG output
  ###
  par(mfrow = c(3,2))
  #imba by logratio from hscrs
  plot1v2(a1 = 2*abs(ab$A1/(ab$A1+ab$A2) - 0.5), a2=log2(ab$A1+ab$A2), mat = ab, clust = clust, cl = cl, w = w,
          xlab = 'log2(Total Copy Number)', ylab = 'Alleleic imbalance 2*|BAF - 0.5|')
  
  #A1 vs A2
  plot1v2(a1 = ab$A1, a2=ab$A2, mat = ab, clust = clust, cl = cl, w = w)
  
  ### Plots of histograms
  myab = histCol(mat = ab, cl = cl, clust = clust, xlab = 'Relative CN homologues pooled')
  myab = histCol(toplot = 'A1', cl = cl, clust=clust, myab = myab, xlab = 'Relative CN homologue 1')
  frame()
  myab = histCol(toplot = 'A2', cl = cl, clust=clust, myab = myab, xlab = 'Relative CN homologue 2')
  
   
  title(main = info$sample.name[i], outer = T, line = -1)
  dev.off()
}
#dev.off()



################################################################
################################################################
#
# Transform
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
library(Rsolnp)

i = 23
dhapseg = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
dhapseg = paste(dhapseg,'/', info$sample.name[i], sep = '')
#load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
seg.dat.fn = file.path(dhapseg, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = ''))
seg.dat = HAPSEGMakeSegObj(seg.dat.fn, gender = 'Female', filter_segs=FALSE, min_probes=10, max_sd=100, 
                           verbose=TRUE)  
as.seg.dat = seg.dat$as.seg.dat
nr = as.integer(nrow(as.seg.dat))
mat = as.seg.dat[1:(nr/2),]
names(mat)[6] = 'A1'
mat$A2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
mat$W = mat$W*2
  
###Read selected segments
tmp = getKnownClustersAndOverlaps(mat = mat, d = file.path(getwd(),'ABSOLUTE', 'SampleData', info$sample.name[i]))
cl = tmp$cl
clust = tmp$clust
w = tmp$w

#Define data for optimization
y = mat$A1
x = mat$A2
xclust = rep(NA, nrow(mat))
yclust = rep(NA, nrow(mat))
cv = read.delim(file.path(getwd(), 'ABSOLUTE', 'SampleData', 'ClusterCNs.txt'))
cv = cv[cv$sample == i,]

ix = which(clust == 'cluster hom')
if (length(ix)>0) xclust[ix] = cv$xhom
if (length(ix)>0) yclust[ix] = cv$yhom

ix = which(clust == 'cluster het')
if (length(ix)>0) xclust[ix] = cv$xhet
if (length(ix)>0) yclust[ix] = cv$yhet

ix = which(clust == 'cluster 1')
if (length(ix)>0) xclust[ix] = cv$x1
if (length(ix)>0) yclust[ix] = cv$y1

ix = which(clust == 'cluster 2')
if (length(ix)>0) xclust[ix] = cv$x2
if (length(ix)>0) yclust[ix] = cv$y2

ix = which(clust == 'cluster 3')
if (length(ix)>0) xclust[ix] = cv$x3
if (length(ix)>0) yclust[ix] = cv$y3

ix = which(clust == 'cluster 4')
if (length(ix)>0) xclust[ix] = cv$x4
if (length(ix)>0) yclust[ix] = cv$y4


#Save the cluster data to use
dd = data.frame(x=x, y=y, xclust = xclust, yclust=yclust, xw=w/sum(w), yw = w/sum(w),
                myclust = clust)
dd = subset(dd, !(is.na(xclust) & is.na(yclust)))


##########################################################
#Optimize one issue at the time


### Scale x
###
dx = subset(dd, !is.na(xclust))
pars = solnp(pars = c(1,0), fun = minx, myclust = dx$myclust,
             X=dx$x, Y=dx$y, xclust = dx$xclust, xw = dx$xw, 
             eqfun = eqfunx, eqB = mean(dx$x), LB = c(0, -1000))$pars
a = pars[1]
b = pars[2]
mat$A2new = a*mat$A2 + b*mat$A1
dd$x = a*dd$x + b*dd$y
param = c(a, b)

### Scale y
###
dy = subset(dd, !is.na(yclust))
pars = solnp(pars = c(0,1), fun = miny, myclust = dy$myclust,
             X=dy$x, Y=dy$y, yclust = dy$yclust, yw = dy$yw, 
             eqfun = eqfuny, eqB = mean(dy$y), LB = c(-1000, 0))$pars
c = pars[1]
d = pars[2]
mat$A1new = c*mat$A2new + d*mat$A1
dd$y = c*dd$x + d*dd$y
param = c(param, c, d)

### Scale  so that x and y fit together
###

pars = solnp(pars = c(1, 1), fun = minxy, myclust = dd$myclust, mat = mat,
             eqfun = eqfuncxy, eqB = 1,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
B = pars[1]
D = pars[2]
mat$A2last =  B*mat$A2new
mat$A1last =  D*mat$A1new
param = c(param, B, D)
names(param) = c('a','b','c','d', 'B','D')
dd$x = B*dd$x
dd$y = D*dd$y


###########################################################
# Figures of outcome
plot1v2(a1 = mat$A1, a2 = mat$A2, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
d = file.path(getwd(), 'ABSOLUTE','SampleData',info$sample.name[i], 'Transform')
dir.create(d, recursive = T, showWarnings = F)
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2orig.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledX.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1new, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledXY.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1last, a2 = mat$A2last, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledAll.png', sep = '')), 
          width=640, height=640)



#Histograms
d = file.path(getwd(), 'ABSOLUTE','SampleData',info$sample.name[i], 'Transform')
myab = histCol(mat = mat, cl = cl, clust = clust, xlab = 'Relative CN (homologues pooled)')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Historig.png', sep = '')), 
          width=640, height=640)

myab = histCol(toplot = c('A1last','A2last'), myab = myab, cl = cl, xlab = 'Relative CN (homologues pooled)')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_HistscaledAll.png', sep = '')), 
          width=640, height=640)

###########################################################


  

#Reset attenuation constant
seg.dat$error.model$AT = 0

#New version or hscr values in as.seg.dat
as.seg.dat$which = 1
as.seg.dat$which[(nr/2+1):nr] = 2
upper = as.seg.dat[as.seg.dat$which == 1,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last')]
upper = merge(upper, tmp, all = T)
upper$copy_num = upper$A1last
upper = upper[order(upper$seg.ix),1:9]
lower = as.seg.dat[as.seg.dat$which == 2,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A2last')]
lower = merge(lower, tmp, all = T)
lower$copy_num = lower$A2last
lower = lower[order(lower$seg.ix),1:9]
as.seg.dat = rbind(upper, lower)
seg.dat$as.seg.dat = as.seg.dat


#New version of allele.segs
allele.segs = seg.dat$allele.segs
allele.segs = as.data.frame(allele.segs)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
allele.segs = merge(allele.segs, tmp, sort = F)
allele.segs$A1.Seg.CN = allele.segs$A1last
allele.segs$A2.Seg.CN = allele.segs$A2last
allele.segs = as.matrix(allele.segs[,1:9])
seg.dat$allele.segs = allele.segs

#New version of seg.info
seg.info = seg.dat$seg.info
seg.info = as.data.frame(seg.info)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
seg.info = merge(seg.info, tmp, sort = F)
seg.info$cn = (seg.info$A1last + seg.info$A2last)/2
seg.info$copy_num = seg.info$cn
seg.info = as.matrix(seg.info[,1:8])
seg.dat$seg.info = seg.info

#Now switch over so that A2 is always the largest (keep index of switched ones)
ix = which(seg.dat$allele.segs[,'A1.Seg.CN'] > seg.dat$allele.segs[,'A2.Seg.CN'])
a1 = seg.dat$allele.segs[,'A1.Seg.CN']
a2 = seg.dat$allele.segs[,'A2.Seg.CN']
seg.dat$allele.segs[,'A1.Seg.CN'] = pmin(a1, a2)
seg.dat$allele.segs[,'A2.Seg.CN'] = pmax(a1, a2)
seg.dat$switched.allele.segs = ix

nr = nrow(seg.dat$as.seg.dat)
a1 = seg.dat$as.seg.dat[1:(nr/2),'copy_num']
a2 = seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num']
seg.dat$as.seg.dat[1:(nr/2),'copy_num'] = pmin(a1, a2)
seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num'] = pmax(a1, a2)

#Print the new object
save(seg.dat, file = file.path(dhapseg, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled_segdat.RData', sep = '')))



################################################################
################################################################
#
# Run ABSOLUTE with scaled haplotype copy values
#
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')



DoAbsolute <- function(myarray, mysample, i, h) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG', h, '/', mysample,'/',mysample,'_',myarray,
                      "_scaled_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute2(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 2, max.ploidy = 15, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample, b.res = 0.125, d.res = 0.125)
  sink()
}

hapsegind = c(2, 1, 1, 1, 1, 2)
names(hapsegind)=as.character(c(21, 27, 40, 51, 20, 23))
  
  foreach (i = c(21, 27, 40, 51, 20, 23), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i, h = hapsegind[as.character(i)])
}



d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute', sep = '')
obj.name = '15,20,21,23,27,40,51'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                     info$sample.name[c(15,20,21,23,27,40,51)], '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)

obj.name = '15,21'
d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute', sep = '')
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "IL", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "AfterScaling", copy_num_type="allelic")


obj.name = '15_subset'
myCreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', info$sample.name[c(15)], 
                                                             '.ABSOLUTE.RData', sep = ''), 
                     indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                     plot.modes = TRUE, verbose=TRUE, myModes = c(1, 2, 4, 10))



#######Junk
foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}

tmp = seg.dat[['obs.scna']][['segtab']]
tmp = tmp[tmp$Chromosome == 3,]
plot(range(c(tmp$Start.pb, tmp$End.bp)), range(tmp$copy_num), type = 'n', 
     xlim = c( 178910000,  178950000))
for (r in 1:nrow(tmp)) lines(c(tmp$Start.bp[r], tmp$End.bp[r]), rep(tmp$copy_num[r], 2))
abline(v = 178927976, col = 'red')

plot(range(c(tmp$tcnStart, tmp$tcnEnd)), range(tmp$tcnMean, na.rm = T), type = 'n', 
     xlim = c( 178910000,  178950000))
for (r in 1:nrow(tmp)) lines(c(tmp$tcnStart[r], tmp$tcnEnd[r]), rep(tmp$tcnMean[r], 2))
abline(v = 178927976, col = 'red')


################################################################
################################################################
#
# BAF evaluation
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
library(IRanges)
library(Rsolnp)
library(limma)

i = 15
solutions = c(1,2,4,10)
d = file.path(getwd(), 'ABSOLUTE', 'ScaledAbsolute')
d = paste(d,'/', '15,21', sep = '')
load(file.path(d, "15,21.PP-modes.data.RData"))
tab = segobj.list[[info$sample.name[i]]][['mode.res']][['mode.tab']][solutions,]
tab = tab[order(tab[,'alpha']),]
as.seg.dat = segobj.list[[info$sample.name[i]]][['as.seg.dat']]
nr = as.integer(nrow(as.seg.dat))
mat = as.seg.dat[1:(nr/2),]
names(mat)[6] = 'A1'
mat$A2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
mat$W = mat$W*2
mat$baf = mat$A2/(mat$A1 + mat$A2)

###Read selected segments
tmp = getKnownClustersAndOverlaps(mat = mat, d = file.path(getwd(),'ABSOLUTE', 'SampleData', info$sample.name[i]))
cl = tmp$cl
clust = tmp$clust
w = tmp$w

#Lowest CN guess of chosen segments
y = mat$A1
x = mat$A2
xclust = rep(NA, nrow(mat))
yclust = rep(NA, nrow(mat))
cv = read.delim(file.path(getwd(), 'ABSOLUTE', 'SampleData', 'ClusterCNs.txt'))
cv = cv[cv$sample == i,]

ix = which(clust == 'cluster hom')
if (length(ix)>0) xclust[ix] = cv$xhom
if (length(ix)>0) yclust[ix] = cv$yhom

ix = which(clust == 'cluster het')
if (length(ix)>0) xclust[ix] = cv$xhet
if (length(ix)>0) yclust[ix] = cv$yhet

ix = which(clust == 'cluster 1')
if (length(ix)>0) xclust[ix] = cv$x1
if (length(ix)>0) yclust[ix] = cv$y1

ix = which(clust == 'cluster 2')
if (length(ix)>0) xclust[ix] = cv$x2
if (length(ix)>0) yclust[ix] = cv$y2

ix = which(clust == 'cluster 3')
if (length(ix)>0) xclust[ix] = cv$x3
if (length(ix)>0) yclust[ix] = cv$y3

ix = which(clust == 'cluster 4')
if (length(ix)>0) xclust[ix] = cv$x4
if (length(ix)>0) yclust[ix] = cv$y4


#Save lowest guesses in dataframe
dd = data.frame(x=x, y=y, xclust = xclust, yclust=yclust, xw=w/sum(w), yw = w/sum(w),
                myclust = clust)
#dd = subset(dd, !(is.na(xclust) & is.na(yclust)))

##########################################################
#
#Segments to use
ix = which(clust != '' & dd$xclust != dd$yclust)

#Observed BAFs
BAFo = mat$baf[ix]

#Weights
W = w[ix]

ss = numeric(nrow(tab))
for (j in 1:nrow(tab)){
  
  #Theoretical BAFs
  add = j-1
  a2 = dd$xclust[ix] + add
  a1 = dd$yclust[ix] + add
  alpha = tab[j,'alpha']
  BAFe = (alpha*a2 + (1-alpha))/(alpha*(a1+a2) + (1-alpha)*2)
  
  #Out
  ss[j] = sum(W*(BAFo-BAFe)**2)
}

plot(tab[,'alpha'], ss, type= 'b', log = 'y')


BAFe = (alpha*a2 + (1-alpha))/(alpha*(a1+a2) + (1-alpha)*2)

ebaf = function(a1, a2, alpha){(alpha*a2 + (1-alpha))/(alpha*(a1+a2) + (1-alpha)*2)}

ebaf(1,2, .28)
ebaf(2,3, .38)


##########################################################
#Optimize one issue at the time


### Scale x
###
dx = subset(dd, !is.na(xclust))
pars = solnp(pars = c(1,0), fun = minx, myclust = dx$myclust,
             X=dx$x, Y=dx$y, xclust = dx$xclust, xw = dx$xw, 
             eqfun = eqfunx, eqB = mean(dx$x), LB = c(0, -1000))$pars
a = pars[1]
b = pars[2]
mat$A2new = a*mat$A2 + b*mat$A1
dd$x = a*dd$x + b*dd$y
param = c(a, b)

### Scale y
###
dy = subset(dd, !is.na(yclust))
pars = solnp(pars = c(0,1), fun = miny, myclust = dy$myclust,
             X=dy$x, Y=dy$y, yclust = dy$yclust, yw = dy$yw, 
             eqfun = eqfuny, eqB = mean(dy$y), LB = c(-1000, 0))$pars
c = pars[1]
d = pars[2]
mat$A1new = c*mat$A2new + d*mat$A1
dd$y = c*dd$x + d*dd$y
param = c(param, c, d)

### Scale  so that x and y fit together
###

pars = solnp(pars = c(1, 1), fun = minxy, myclust = dd$myclust, mat = mat,
             eqfun = eqfuncxy, eqB = 1,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
B = pars[1]
D = pars[2]
mat$A2last =  B*mat$A2new
mat$A1last =  D*mat$A1new
param = c(param, B, D)
names(param) = c('a','b','c','d', 'B','D')
dd$x = B*dd$x
dd$y = D*dd$y


###########################################################
# Figures of outcome
plot1v2(a1 = mat$A1, a2 = mat$A2, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
d = file.path(getwd(), 'ABSOLUTE','SampleData',info$sample.name[i], 'Transform')
dir.create(d, recursive = T, showWarnings = F)
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2orig.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledX.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1new, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledXY.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1last, a2 = mat$A2last, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledAll.png', sep = '')), 
          width=640, height=640)



#Histograms
d = file.path(getwd(), 'ABSOLUTE','SampleData',info$sample.name[i], 'Transform')
myab = histCol(mat = mat, cl = cl, clust = clust, xlab = 'Relative CN (homologues pooled)')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Historig.png', sep = '')), 
          width=640, height=640)

myab = histCol(toplot = c('A1last','A2last'), myab = myab, cl = cl, xlab = 'Relative CN (homologues pooled)')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_HistscaledAll.png', sep = '')), 
          width=640, height=640)

###########################################################




#Reset attenuation constant
seg.dat$error.model$AT = 0

#New version or hscr values in as.seg.dat
as.seg.dat$which = 1
as.seg.dat$which[(nr/2+1):nr] = 2
upper = as.seg.dat[as.seg.dat$which == 1,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last')]
upper = merge(upper, tmp, all = T)
upper$copy_num = upper$A1last
upper = upper[order(upper$seg.ix),1:9]
lower = as.seg.dat[as.seg.dat$which == 2,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A2last')]
lower = merge(lower, tmp, all = T)
lower$copy_num = lower$A2last
lower = lower[order(lower$seg.ix),1:9]
as.seg.dat = rbind(upper, lower)
seg.dat$as.seg.dat = as.seg.dat


#New version of allele.segs
allele.segs = seg.dat$allele.segs
allele.segs = as.data.frame(allele.segs)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
allele.segs = merge(allele.segs, tmp, sort = F)
allele.segs$A1.Seg.CN = allele.segs$A1last
allele.segs$A2.Seg.CN = allele.segs$A2last
allele.segs = as.matrix(allele.segs[,1:9])
seg.dat$allele.segs = allele.segs

#New version of seg.info
seg.info = seg.dat$seg.info
seg.info = as.data.frame(seg.info)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
seg.info = merge(seg.info, tmp, sort = F)
seg.info$cn = (seg.info$A1last + seg.info$A2last)/2
seg.info$copy_num = seg.info$cn
seg.info = as.matrix(seg.info[,1:8])
seg.dat$seg.info = seg.info

#Now switch over so that A2 is always the largest (keep index of switched ones)
ix = which(seg.dat$allele.segs[,'A1.Seg.CN'] > seg.dat$allele.segs[,'A2.Seg.CN'])
a1 = seg.dat$allele.segs[,'A1.Seg.CN']
a2 = seg.dat$allele.segs[,'A2.Seg.CN']
seg.dat$allele.segs[,'A1.Seg.CN'] = pmin(a1, a2)
seg.dat$allele.segs[,'A2.Seg.CN'] = pmax(a1, a2)
seg.dat$switched.allele.segs = ix

nr = nrow(seg.dat$as.seg.dat)
a1 = seg.dat$as.seg.dat[1:(nr/2),'copy_num']
a2 = seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num']
seg.dat$as.seg.dat[1:(nr/2),'copy_num'] = pmin(a1, a2)
seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num'] = pmax(a1, a2)

#Print the new object
save(seg.dat, file = file.path(dhapseg, paste(info$sample.name[i], '_', info$arrname[i], 
                                              '_scaled_segdat.RData', sep = '')))





################################################################
################################################################
#
# Simulation of dataset for Hapseg
# This never passed HAPSEG
################################################################
################################################################


info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample = substr(info$Sample.ID, 1, 4)
info$sample.name = paste(rownames(info), info$sample, sep = '.')
info$myarray = info$arrname


i = 15

#Read real data and segments
d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i])
data = read.delim(file.path(d, 'Total,fracB.txt'), stringsAsFactors = F)
d = file.path(getwd(), 'ABSOLUTE', 'SampleData', info$sample.name[i], 'Segdat.txt')
newd = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'SampleData', info$sample.name[i])
dir.create(newd, recursive = T, showWarnings = F)
file.copy(d, file.path(newd, 'Segdat.txt'))
seg = read.delim(d)

#Build matrix of possible copy numbers
cns = rbind(data.frame(t = 6, maj = 3:4, min = 3:2),
            data.frame(t = 7, maj = 4:5, min = 3:2),
            data.frame(t = 8, maj = 4:6, min = 4:2),
            data.frame(t = 9, maj = 5, min = 4))

#Simulate observations from those copy numbers
set.seed(100)
copyrow = ceiling(runif(nrow(seg))*nrow(cns))
cop = cns[copyrow,]
cop$baf1 = cop$min/cop$t
cop$baf2 = cop$maj/cop$t
cop$imba = cop$baf2 - cop$baf1

data$total = NA
data$fracB = NA
data$A = NA
data$B = NA
for (r in 1:nrow(seg)){
  ix = which(data$chromosome == seg$Chromosome[r] & data$position>=seg$Start.Position[r] & 
               data$position<=seg$End.Position[r])
  noise = matrix(rnorm(2*length(ix), sd = .5), ncol = 2)
  total = (cop$t[r] + noise[,1] + noise[,2])
  baf1 = (noise[,1])/total
  baf2 = (cop$min[r] + noise[,1])/total
  baf3 = (cop$maj[r] + noise[,1])/total
  baf4 = (cop$t[r] + noise[,1])/total
  rand = ceiling(runif(length(ix)))
  A = numeric(length(ix))
  B = numeric(length(ix))
  fracB = numeric(length(ix))
  for (j in 1:4){
    ind = which(rand == j)
    bafj = get(paste('baf', j, sep = ''))
    A[ind] = total[ind]*(1-bafj[ind])
    B[ind] = total[ind]*bafj[ind]
    fracB[ind] = bafj[ind]
  }
  data$A[ix] = A
  data$B[ix] = B
  data$fracB[ix] = fracB
  data$total[ix] = total
}



#Example BAF plot
plot(data$position[ix], baf1, pch = '.', ylim = c(0,1), xlab =paste( 'Position on Chromosome',
          seg$Chromosome[r]), ylab = 'B Allele Fraction')
points(data$position[ix], baf2, pch = '.')
points(data$position[ix], baf3, pch = '.')
points(data$position[ix], baf4, pch = '.')
dev.print(png, file = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'SampleData', 
         info$sample.name[i], 'BAFs.png'), height = 400, width = 640)


#Print dataset
data$baf = pmax(data$fracB, 0)
data$baf = pmin(data$fracB, 1)
data$A = data$total*(1-data$baf)
data$B = data$total*(data$baf)

dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'SampleData', info$sample.name[i])
dir.create(d, recursive = T, showWarnings = F)
save(dat, file = file.path(d, 'snpfile.RData'))

dat = subset(data, select = c('unitNames','chromosome','position','total', 'fracB'))
write.table(data, file.path(d, 'Total,fracB.txt'), row.names = F, quote = F, sep = '\t')



################################################################
################################################################
#
# Sample 15 rescale again
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
library(IRanges)
library(Rsolnp)
library(limma)

i = 21
d = file.path(getwd(), 'ABSOLUTE', 'ScaledAbsolute')
d = paste(d,'/', '15,20,21,23,27,40,51', sep = '')
load(file.path(d, "15,20,21,23,27,40,51.PP-modes.data.RData"))
modpars = segobj.list[['21.4669']]$mode.res$mode.tab[2,]
delta = modpars['delta']

d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
d = paste(d,'/', info$sample.name[i], sep = '')
seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_scaled_segdat.RData', sep = ''))
seg.dat = myMakeSegObj(seg.dat.fn, filter_segs=FALSE, min_probes=10, max_sd=100, 
                       verbose=TRUE)  
as.seg.dat = seg.dat$as.seg.dat



#############################################################
#Rescale and print new dataset

as.seg.dat$copy_num = as.seg.dat$copy_num + 4*delta
e.cr = sum(as.seg.dat$copy_num* as.seg.dat$W)
as.seg.dat$copy_num = as.seg.dat$copy_num/e.cr
e.cr = sum(as.seg.dat$copy_num* as.seg.dat$W)
seg.dat$as.seg.dat = as.seg.dat



#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', info$sample.name[i])
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled2_segdat.RData', sep = '')))

################################################################
################################################################
#
# Run ABSOLUTE with scaled2 haplotype copy values (sample 15)
#
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')



DoAbsolute <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                      "_scaled2_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/Scaled2Absolute', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/Scaled2Absolute/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute2(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 2, max.ploidy = 12, 
                 primary.disease = "Breast Cancer", 
                 platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                 max.non.clonal = 0, 
                 max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                 min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}

foreach (i = c(15), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}



d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute', sep = '')
obj.name = '15,21'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                                                           info$sample.name[c(15,21)], '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)

obj.name = '15,21'
d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute', sep = '')
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "IL", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "AfterScaling", copy_num_type="allelic")


obj.name = '15_subset'
myCreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', info$sample.name[c(15)], 
                                                             '.ABSOLUTE.RData', sep = ''), 
                     indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                     plot.modes = TRUE, verbose=TRUE, myModes = c(1:8,12))





################################################################
################################################################
#
# Simulation of data and Lilekihood study with ABSOLUTE
# 
################################################################
################################################################


library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')


i = 21
myarray = info$myarray[i] 
mysample = info$sample.name[i]
sample.name <- mysample
maf.fn = NULL
seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                    "_scaled_segdat.RData", sep="")
results.dir <- paste(getwd(), '/ABSOLUTE/Simulation/ScaledAbsolute', sep = '')
print(paste("Starting myarray", mysample, "at", results.dir))
log.dir <- paste(getwd(), '/ABSOLUTE/Simulation/ScaledAbsolute/abs_logs', sep = '') 
if (!file.exists(log.dir)) {
  dir.create(log.dir, recursive=TRUE)
}
if (!file.exists(results.dir)) {
  dir.create(results.dir, recursive=TRUE)
}
sigma.p = 0
max.sigma.h = 0.01
min.ploidy = 2
max.ploidy = 12 
primary.disease = "Breast Cancer"
platform = 'SNP_6.0'
max.as.seg.count = 1500
max.non.clonal = 0
max.neg.genome = 0
copy_num_type = 'allelic'
verbose=TRUE
maf.fn = maf.fn
min.mut.af = 0.03
output.fn.base = mysample
b.res = 0.125
d.res = 0.125

#MyRunabsolute2
kQ = 15
kTauDom = c(min.ploidy, max.ploidy)
kSigmaHDom = c(0, max.sigma.h)
kPiSomThetaQ = c(100, 50, rep(2, (kQ - 2)))
mut_class_w = list(SM = 0.5, GL = 0, SC = 0.5, OL = 0.001, 
                   Pi_SM = 15, Pi_SC = 15)
platform = match.arg(platform, c("SNP_6.0", "Illumina_WES", 
                                 "SNP_250K_STY"))
if (platform %in% c("SNP_6.0", "SNP_250K_STY")) {
  filter_segs = FALSE
} else if (platform == "Illumina_WES") {
  filter_segs = TRUE
} else {
  stop("Unsupported platform: ", platform)
}
min_probes = 10
max_sd = 100
copy_num_type = match.arg(copy_num_type, c("allelic", "total"))
if (copy_num_type == "total") {
  pi_theta_qz = list(W = c(1, 25, 100, 25, 10, 5, rep(1, 
                                                      kQ - 6), 10), PC = 0.05)
  ABSOLUTE:::set_total_funcs()
} else if (copy_num_type == "allelic") {
  pi_theta_qz = list(W = c(25, 100, 25, 10, 5, rep(1, kQ - 
                                                     5), 10), PC = 0.05)
  ABSOLUTE:::set_allelic_funcs()
}  else {
  stop("Unsupported copy number type: ", copy_num_type)
}
SubclonalPost <<- ABSOLUTE:::UnifSubclonalPost
sigma.h = sigma.p
tmp.dir = file.path(results.dir, "tmp")
dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
seg.dat = myMakeSegObj(seg.dat.fn, min_probes = min_probes, 
                       max_sd = max_sd, filter_segs = filter_segs, verbose = verbose)
#Simulate new a1 and a2 levels

#For sample 15
mydelta = 0.12873490
myb = seq(0, by = 2*mydelta, length.out = 3)
hist(seg.dat$as.seg.dat$copy_num, breaks = 200)
mu = seq(.75, by = .25, length.out = 3)
sig = 0.01
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample15', 'Simulation15.RData'))
#  obs[hetindex == 1] = runif(table(hetindex)[2])*3 + .75 #high
#  obs[hetindex == 1] = runif(table(hetindex)[2]) #low
#  obs[hetindex == 1] = runif(table(hetindex)[2])/2 + .75 #middle

#For sample 21 version 1
mydelta = 0.12493021
#mydelta = 0.12745397
myb = seq(0.6, by = -2*mydelta, length.out = 3)
hist(seg.dat$as.seg.dat$copy_num, breaks = 200)
mu = seq(.6, by = .25, length.out = 3)
sig = 0.05
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample21', 'Simulation21.RData'))
#mode.tabspar = mode.res$mode.tab

#For sample 21 version 2
#mydelta = 0.12493021
mydelta = 0.12745397
myb = seq(0.6, by = -2*mydelta, length.out = 3)
#hist(seg.dat$as.seg.dat$copy_num, breaks = 200)
mu = seq(.6, by = .25, length.out = 3)
sig = 0.05
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample21', 'Simulation21_2.RData'))
#mode.tabspar = mode.res$mode.tab

#load("/wehisan/home/allstaff/l/lonnstedt/RESPONSIFY/Simulation.RData")
#load("/wehisan/home/allstaff/l/lonnstedt/RESPONSIFY/ABSOLUTE/Simulation/ScaledAbsolute/Sample21/Simulation21.RData")


LLs = matrix(NA, ncol = 3, nrow = 50)
for (sim in 1:50){
  means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
  obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
  hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
#  obs[hetindex == 1] = runif(table(hetindex)[2])*3 + .6 #high
#  obs[hetindex == 1] = runif(table(hetindex)[2]) #low
#  obs[hetindex == 1] = runif(table(hetindex)[2])/2 + .6 #middle
  obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
  #hist(obs, breaks = 100, xlim = c(0, 2))
  seg.dat$as.seg.dat$copy_num = obs
  
#  hist(seg.dat$as.seg.dat$copy_num, breaks = 1000, xlim = c(0, 2))
#  hist(obs, breaks = 1000, xlim = c(0, 2))
  
  
  #Run absolute
#  seg.dat[["primary.disease"]] = primary.disease
#  seg.dat[["group"]] = ABSOLUTE:::DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  seg = seg.dat[["obs.scna"]][["segtab"]]
  e.cr = sum(seg[, "W"] * seg[, "copy_num"])
  if (verbose) {
    print(paste("Expected copy-ratio = ", round(e.cr, 5), 
                sep = ""))
  }
  mode.res = list(mode.flag = NA)
  if (nrow(seg) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  if ((e.cr < 0.75) || (e.cr > 1.25)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
  ##########
  ##########
  maf = NULL
#  data(ChrArmsDat, package = "ABSOLUTE")
    mut = NA
  #myFitSample
  seg.obj = seg.dat
  Q = kQ
  tau.dom = kTauDom
  sigma.h.dom = kSigmaHDom
  
  kThetaQz = rep(1/(Q + 1), Q + 1)
  obs = seg.obj[["obs.scna"]]
  
  b = myb[3:1]
  delta = rep(mydelta, length(b))
  alpha = 2*delta[1]/(2*delta[1] + b)
  tau = (1-b)/delta[1]
  AT = rep(0, length(b))
  LL = rep(NA, length(b))
  mode_curv = rep(NA, length(b))
  mode.tab = cbind(alpha, tau, AT, b, delta, LL, mode_curv)
  
  if (is.na(mode.tab)) {
    return(list(mode.flag = "ERROR"))
  }
  n.modes = nrow(mode.tab)
  theta.qz.hat = matrix(NA, nrow = n.modes, ncol = Q + 1)
  theta.q.tab = array(NA, dim = c(n.modes, Q))
  if (verbose) {
    print(paste("Optimizing LL(data, theta.qz, sigma.h | comb) for ", 
                n.modes, " modes: ", sep = ""))
  }
  seg.z.tab = array(NA, dim = c(n.modes, length(obs[["W"]])))
  seg.qz.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                 Q + 1))
  seg.q.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                Q))
  if (obs[["data.type"]] == "ALLELIC") {
    ab.tab = array(NA, dim = c(n.modes, Q + 1))
    chr.arm.tab = array(NA, dim = c(n.modes, 2, nrow(chr.arms.dat), 
                                    Q))
  }
  if (obs[["data.type"]] == "TOTAL") {
    ab.tab = NULL
    chr.arm.tab = array(NA, dim = c(n.modes, 1, nrow(chr.arms.dat), 
                                    Q))
  }
  dimnames(chr.arm.tab)[[3]] = rownames(chr.arms.dat)
  mode.tab = cbind(mode.tab, NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mode.tab)[c((ncol(mode.tab) - 7):ncol(mode.tab))] = c("genome mass", 
                                                                 "sigma.h.hat", "theta.z.hat", "frac.het", "log.partition", 
                                                                 "entropy", "clust_LL", "post_LL")
  mode.tab = cbind(mode.tab, somatic.mut.ll = NA)
  for (i in 1:n.modes) {
    delta = mode.tab[i, "delta"]
    b = mode.tab[i, "b"]
    obs[["error.model"]][["fit.at"]] = mode.tab[i, "AT"]
    comb = GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
    res = ABSOLUTE:::OptThetaQzSigmaH(obs, comb, sigma.h.dom, pi_theta_qz, 
                                      verbose = verbose)
    mode.tab[i, "sigma.h.hat"] = res[["sigma.h.hat"]]
    mode.tab[i, "log.partition"] = res[["LL"]]
    mode.tab[i, "post_LL"] = res[["LL"]]
  }
  LLs[sim,] = mode.tab[,'log.partition']
}
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation',  'ScaledAbsolute', 'Sample21', 'Highhet21_2.RData'))
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation',  'ScaledAbsolute', 'Sample21', 'Lowhet21_2.RData'))
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation',  'ScaledAbsolute', 'Sample21', 'Middlehet21_2.RData'))
save.image(file.path(getwd(), 'ABSOLUTE', 'Simulation',  'ScaledAbsolute', 'Sample21', 'Nohet21_2.RData'))


################
# No heterogeneity

d = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample21')
namn = '21_2'
lab = c('0, 1, 2', '1, 2, 3','2, 3, 4')
#lab = c('1, 2, 3','2, 3, 4', '3, 4, 5')

load(file.path(d, paste('Nohet',namn,'.RData', sep = '')))


#Hist
m = apply(LLs, 1, which.max)
hist(m, xlim = c(0.5,3.5), col = 'salmon', xlab = 'Solution', main= 'ML solution', xaxt = 'n', breaks = 10)
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_hist_nohet', namn, '.png', sep = '')), width = 340, height = 340)

#Lines
matplot(LLs, type = 'n', xlim = c(0.5,3.5), xaxt = 'n', ylab = 'Log likelihood', xlab = 'Solution')
axis(1, at = 1:3, label = lab)
grid(col = 1)
for (j in 1:nrow(LLs)){
  lines(LLs[j,], col = j)
}
dev.print(png, file = file.path(d, paste('LLs_nohet', namn, '.png', sep = '')), width = 340, height = 540)

#Boxplot
boxplot(LLs, xlim = c(0.5,3.5), col = 'darkseagreen', xlab = 'Solution', main= '', xaxt = 'n',
        ylab = 'Log lokelihood')
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_boxplot_nohet', namn, '.png', sep = '')), width = 340, height = 340)


################
# High CN heterogeneity

load(file.path(d, paste('Highhet',namn,'.RData', sep = '')))
LLs = LLs[1:50,]

#Hist
m = apply(LLs, 1, which.max)
hist(m, xlim = c(0.5,3.5), col = 'salmon', xlab = 'Solution', main= 'ML solution', xaxt = 'n',
     breaks = 10)
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_hist_highhet', namn, '.png', sep = '')), width = 340, height = 340)

#Lines
matplot(LLs, type = 'n', xlim = c(0.5,3.5), xaxt = 'n', ylab = 'Log likelihood', xlab = 'Solution')
axis(1, at = 1:3, label = lab)
grid(col = 1)
for (j in 1:nrow(LLs)){
  lines(LLs[j,], col = j)
}
dev.print(png, file = file.path(d, paste('LLs_highhet', namn, '.png', sep = '')), width = 340, height = 540)

#Boxplot
boxplot(LLs, xlim = c(0.5,3.5), col = 'darkseagreen', xlab = 'Solution', main= '', xaxt = 'n',
        ylab = 'Log lokelihood')
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_boxplot_highhet', namn, '.png', sep = '')), width = 340, height = 340)


################
# Low CN heterogeneity

load(file.path(d, paste('Lowhet',namn,'.RData', sep = '')))
LLs = LLs[1:50,]

#Hist
m = apply(LLs, 1, which.max)
hist(m, xlim = c(0.5,3.5), col = 'salmon', xlab = 'Solution', main= 'ML solution', xaxt = 'n',
     breaks = 10)
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_hist_lowhet', namn, '.png', sep = '')), width = 340, height = 340)

#Lines
matplot(LLs, type = 'n', xlim = c(0.5,3.5), xaxt = 'n', ylab = 'Log likelihood', xlab = 'Solution')
axis(1, at = 1:3, label = lab)
grid(col = 1)
for (j in 1:nrow(LLs)){
  lines(LLs[j,], col = j)
}
dev.print(png, file = file.path(d, paste('LLs_lowhet', namn, '.png', sep = '')), width = 340, height = 540)

#Boxplot
boxplot(LLs, xlim = c(0.5,3.5), col = 'darkseagreen', xlab = 'Solution', main= '', xaxt = 'n',
        ylab = 'Log lokelihood')
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_boxplot_lowhet', namn, '.png', sep = '')), width = 340, height = 340)



################
# Middle CN heterogeneity

load(file.path(d, paste('Middlehet',namn,'.RData', sep = '')))
LLs = LLs[1:50,]

#Hist
m = apply(LLs, 1, which.max)
hist(m, xlim = c(0.5,3.5), col = 'salmon', xlab = 'Solution', main= 'ML solution', xaxt = 'n',
     breaks = 10)
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_hist_middlehet', namn, '.png', sep = '')), width = 340, height = 340)

#Lines
matplot(LLs, type = 'n', xlim = c(0.5,3.5), xaxt = 'n', ylab = 'Log likelihood', xlab = 'Solution')
axis(1, at = 1:3, label = lab)
grid(col = 1)
for (j in 1:nrow(LLs)){
  lines(LLs[j,], col = j)
}
dev.print(png, file = file.path(d, paste('LLs_middlehet', namn, '.png', sep = '')), width = 340, height = 540)

#Boxplot
boxplot(LLs, xlim = c(0.5,3.5), col = 'darkseagreen', xlab = 'Solution', main= '', xaxt = 'n',
        ylab = 'Log lokelihood')
axis(1, at = 1:3, label = lab)
dev.print(png, file = file.path(d, paste('LLs_boxplot_middlehet', namn, '.png', sep = '')), width = 340, height = 340)

################
# Relative CN distributions


#Relative CN histograms sample 15
d = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample15')
load(file.path(d, 'Simulation15.RData'))
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Original data sample 15',
     xlab = 'Relatvie CN', xlim = c(0,3))
dev.print(png, file = file.path(d, 'Relative_CN_original.png'), width = 340, height = 340)


means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_nohet.png'), width = 340, height = 340)


means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .1)
obs[hetindex == 1] = runif(table(hetindex)[2])*3 + .75 #high
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_highhet.png'), width = 340, height = 340)

means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .1)
obs[hetindex == 1] = runif(table(hetindex)[2])/2 + .75 #middle
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_middlehet.png'), width = 340, height = 340)

means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .1)
obs[hetindex == 1] = runif(table(hetindex)[2]) #low
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_lowhet.png'), width = 340, height = 340)



#Relative CN histograms sample 21
d = file.path(getwd(), 'ABSOLUTE', 'Simulation', 'ScaledAbsolute', 'Sample21')
load(file.path(d, 'Simulation21.RData'))
hist(seg.dat$as.seg.dat$copy_num, breaks = 150, main = 'Original data sample 21',
     xlab = 'Relatvie CN', xlim = c(0,3))
dev.print(png, file = file.path(d, 'Relative_CN_original.png'), width = 340, height = 340)


means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_nohet.png'), width = 340, height = 340)


means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
obs[hetindex == 1] = runif(table(hetindex)[2])*3 + .6 #high
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_highhet.png'), width = 340, height = 340)

means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
obs[hetindex == 1] = runif(table(hetindex)[2])/2 + .6 #middle
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_middlehet.png'), width = 340, height = 340)

means = mu[ceiling(runif(nrow(seg.dat$as.seg.dat))*3)]
obs = rnorm(nrow(seg.dat$as.seg.dat), means, sig)
hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
obs[hetindex == 1] = runif(table(hetindex)[2]) #low
obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
seg.dat$as.seg.dat$copy_num = obs
hist(seg.dat$as.seg.dat$copy_num, breaks = 200, main = 'Simulated dataset',
     xlab = 'Relatvie CN', xlim = c(0,3), yaxt = 'n', ylab = '')
dev.print(png, file = file.path(d, 'Relative_CN_lowhet.png'), width = 340, height = 340)





################################################################
################################################################
#
# Simulation 2 of data and Lilekihood study with ABSOLUTE
# 
################################################################
################################################################


library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')


i = 15
myarray = info$myarray[i] 
mysample = info$sample.name[i]
sample.name <- mysample
maf.fn = NULL
seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                    "_scaled_segdat.RData", sep="")
results.dir <- paste(getwd(), '/ABSOLUTE/Simulation/CNbasedAbsolute', sep = '')
print(paste("Starting myarray", mysample, "at", results.dir))
log.dir <- paste(getwd(), '/ABSOLUTE/Simulation/CNbasedAbsolute/abs_logs', sep = '') 
if (!file.exists(log.dir)) {
  dir.create(log.dir, recursive=TRUE)
}
if (!file.exists(results.dir)) {
  dir.create(results.dir, recursive=TRUE)
}
sigma.p = 0
max.sigma.h = 0.01
min.ploidy = 2
max.ploidy = 12 
primary.disease = "Breast Cancer"
platform = 'SNP_6.0'
max.as.seg.count = 1500
max.non.clonal = 0
max.neg.genome = 0
copy_num_type = 'allelic'
verbose=TRUE
maf.fn = maf.fn
min.mut.af = 0.03
b.res = 0.125
d.res = 0.125

#MyRunabsolute2
kQ = 15
kTauDom = c(min.ploidy, max.ploidy)
kSigmaHDom = c(0, max.sigma.h)
kPiSomThetaQ = c(100, 50, rep(2, (kQ - 2)))
mut_class_w = list(SM = 0.5, GL = 0, SC = 0.5, OL = 0.001, 
                   Pi_SM = 15, Pi_SC = 15)
platform = match.arg(platform, c("SNP_6.0", "Illumina_WES", 
                                 "SNP_250K_STY"))
if (platform %in% c("SNP_6.0", "SNP_250K_STY")) {
  filter_segs = FALSE
} else if (platform == "Illumina_WES") {
  filter_segs = TRUE
} else {
  stop("Unsupported platform: ", platform)
}
min_probes = 10
max_sd = 100
copy_num_type = match.arg(copy_num_type, c("allelic", "total"))
if (copy_num_type == "total") {
  pi_theta_qz = list(W = c(1, 25, 100, 25, 10, 5, rep(1, 
                                                      kQ - 6), 10), PC = 0.05)
  ABSOLUTE:::set_total_funcs()
} else if (copy_num_type == "allelic") {
  pi_theta_qz = list(W = c(25, 100, 25, 10, 5, rep(1, kQ - 
                                                     5), 10), PC = 0.05)
  ABSOLUTE:::set_allelic_funcs()
}  else {
  stop("Unsupported copy number type: ", copy_num_type)
}
SubclonalPost <<- ABSOLUTE:::UnifSubclonalPost
sigma.h = sigma.p
tmp.dir = file.path(results.dir, "tmp")
dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
seg.dat = myMakeSegObj(seg.dat.fn, min_probes = min_probes, 
                       max_sd = max_sd, filter_segs = filter_segs, verbose = verbose)
#Simulate new a1 and a2 levels

#load("/wehisan/home/allstaff/l/lonnstedt/RESPONSIFY/Simulation.RData")
#load(file.path(getwd(),"/ABSOLUTE/Simulation/ScaledAbsolute/Sample15/Simulation15.RData"))


for (sim in 1:50){
  output.fn.base = paste(mysample, 'sim', sep = '.')
  alpha = .8
  sig = 0.01
  nr = nrow(seg.dat$as.seg.dat)
  cns = ceiling(runif(nr)*3)+2
  ix1 = 1:(nr/2)
  ix2 = (nr/2+1):nr
  tau = sum(cns[ix1]*seg.dat$as.seg.dat$W[ix1]*2) + sum(cns[ix2]*seg.dat$as.seg.dat$W[ix2]*2)
  err = rnorm(nr, 0, sig)
  CNs = alpha * cns + (1-alpha) * 1
  obs = CNs + err
  
  #hetindex = rbinom(nrow(seg.dat$as.seg.dat), size = 1, prob = .3)
  #  obs[hetindex == 1] = runif(table(hetindex)[2])*3 + .6 #high
  #  obs[hetindex == 1] = runif(table(hetindex)[2]) #low
  #  obs[hetindex == 1] = runif(table(hetindex)[2])/2 + .6 #middle
  obs = obs/sum(seg.dat$as.seg.dat[, "W"] * obs)
  #hist(obs, breaks = 100, xlim = c(0, 2))
  seg.dat$as.seg.dat$copy_num = obs
  
  #  hist(seg.dat$as.seg.dat$copy_num, breaks = 1000, xlim = c(0, 2))
  #  hist(obs, breaks = 1000, xlim = c(0, 2))
  
  
  #Run absolute
  #  seg.dat[["primary.disease"]] = primary.disease
  #  seg.dat[["group"]] = ABSOLUTE:::DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  seg = seg.dat[["obs.scna"]][["segtab"]]
  e.cr = sum(seg[, "W"] * seg[, "copy_num"])
  if (verbose) {
    print(paste("Expected copy-ratio = ", round(e.cr, 5), 
                sep = ""))
  }
  mode.res = list(mode.flag = NA)
  if (nrow(seg) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  if ((e.cr < 0.75) || (e.cr > 1.25)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
   maf = NULL
  #  data(ChrArmsDat, package = "ABSOLUTE")
  mut = NA
  
  
  mode.res = myFitSample(seg.obj = seg.dat, mut, Q = kQ, pi_theta_qz, sigma.h, 
                         tau.dom = kTauDom, sigma.h.dom = kSigmaHDom, 
                         chr.arms.dat, verbose = verbose, d.res = d.res, b.res = b.res)
  mode.res = myapply_subclonal_scna_model(segobj = seg.dat, mode_res = mode.res, 
                                          verbose = verbose)
  if (inherits(mode.res, "try-error")) {
    mode.res = list(mode.flag = "FAIL")
  }
  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = ABSOLUTE:::GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, 
                                        max.non.clonal, max.neg.genome, kQ, verbose = verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag = "ALPHA_TAU_DOM")
    } else {
      mode.res = ABSOLUTE:::ReorderModeRes(mode.res, !bad.ix)
    }
  }
  if (is.na(mode.res[["mode.flag"]])) {
    data(ChrArmPriorDb, package = "ABSOLUTE")
    model.id = ifelse(seg.dat[["group"]] %in% names(train.obj), 
                      seg.dat[["group"]], "Primary")
    mode.res = ABSOLUTE:::ApplyChrArmPrior(mode.res, model.id, train.obj, 
                                           verbose = verbose)
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("maf was not NULL, but min.mut.af was")
      }
      seg.dat[["mut.cn.dat"]] = mut.cn.dat
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, seg.dat$obs_SCNA, 
                                                  mut.cn.dat, kPiSomThetaQ, mut_class_w, kQ, verbose = verbose)
    }
    mode.res = ABSOLUTE:::WeighSampleModes(mode.res)
    mode.res[["call.status"]] = ABSOLUTE:::GetCallStatus(mode.res, seg.dat[["obs.scna"]][["W"]])
  }
  seg.dat[["mode.res"]] = mode.res
  if (is.null(output.fn.base)) {
    output.fn.base = ifelse(is.null(seg.dat$array.name), 
                            sample.name, seg.dat$array.name)
  }
  file.base = paste(output.fn.base, ".ABSOLUTE", sep = "")
  if (is.na(mode.res[["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir, paste(file.base, 
                                                 "plot.pdf", sep = "_"))
    ABSOLUTE:::AbsoluteResultPlot(sample.pdf.fn, seg.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }
  }
  seg.dat$version = 1.1
  save(seg.dat, file = file.path(results.dir, paste(file.base, 
                                                    "RData", sep = ".")))
  
  
  
  
  
  ################################################################
  ################################################################
  #
  # allelic imbalance plots
  # 
  ################################################################
  ################################################################
  ebaf = function(a1, a2, alpha){(alpha*a2 + (1-alpha))/(alpha*(a1+a2) + (1-alpha)*2)}
  etot = function(a1, a2, alpha){(a1 + a2)*alpha + (1-alpha)*2}
  
  par (mfcol = c(1,3), xpd = NA)
  ybar = 1.5
  alphas = c(.28, .39, .64)
  for (i in 1:3){
    alpha = alphas[i]
    a1 = 0:7
    a2 = 0:7
    #alpha = .5
    bafs = outer(a1, a2, alpha, FUN = ebaf)
    ai = 2*abs(bafs-.5)
    dimnames(ai) = list(as.character(a1), as.character(a2))
    tots = outer(a1, a2, alpha, FUN = etot)
    dimnames(tots) = list(as.character(a1), as.character(a2))
    tots[1] = NA
    plot(tots, ai, main  = paste('alpha =', alpha), 
         ylim = c(0,1), xlim = c(0, 10), xlab = 'Total CN',
         ylab = 'Allelic Imbalance')
    points(tots[as.character(i-1),as.character(i+1)], 
           ai[as.character(i-1),as.character(i+1)], col = 'yellow', pch = 16)
    points(tots[as.character(i),as.character(i)], 
           ai[as.character(i),as.character(i)], col = 'pink', pch = 16)
    points(tots[as.character(i),as.character(i+1)], 
           ai[as.character(i),as.character(i+1)], col = 'red', pch = 16)
    points(tots[as.character(i+1),as.character(i+1)], 
           ai[as.character(i+1),as.character(i+1)], col = 'green3', pch = 16)
    points(tots[as.character(i+2),as.character(i+1)], 
           ai[as.character(i+2),as.character(i+1)], col = 'blue', pch = 16)
    abline(h = ai[as.character(i-1),as.character(i+1)], col = 'yellow')
    grid(col = 'grey')
    #totr = tots/(ybar + i -1)
    #plot(totr, ai, main  = paste('alpha =', alpha), 
    #     ylim = c(0,1), xlim = c(0, 4), xlab = 'Relative CN',
    #     ylab = 'Allelic Imbalance')
    #points(totr[as.character(i-1),as.character(i+1)], 
    #       ai[as.character(i-1),as.character(i+1)], col = 'yellow', pch = 16)
    #points(totr[as.character(i),as.character(i)], 
    #       ai[as.character(i),as.character(i)], col = 'pink', pch = 16)
    #points(totr[as.character(i),as.character(i+1)], 
    #       ai[as.character(i),as.character(i+1)], col = 'red', pch = 16)
    #points(totr[as.character(i+1),as.character(i+1)], 
    #       ai[as.character(i+1),as.character(i+1)], col = 'green3', pch = 16)
    #points(totr[as.character(i+2),as.character(i+1)], 
    #       ai[as.character(i+2),as.character(i+1)], col = 'blue', pch = 16)
    #grid(col = 'grey')
  }
  dev.off()
  par(mfcol = c(2,3))
  alphas = c(.28, .39, .64)
  for (i in 1:3){
    alpha = alphas[i]
    a1 = c(0, 1, 1, 2, 2) + i - 1
    a2 = c(2, 1, 2, 2, 3) + i - 1
    a1 = alpha*a1 +(1-alpha)*1
    a2 = alpha*a2 +(1-alpha)*1
    h = hist(c(a1, a2), breaks = 100, xlim = c(0, 3.5), main  = paste('alpha =', alpha),
             xlab = 'Homologue Specific CN')
    abline(v = alpha*c(0:6) + 1-alpha, col = 'green') 
    r = c(a1, a2)/mean(c(a1, a2))
    h = hist(r, breaks = 100, xlim = c(0, 1.5), main  = paste('alpha =', alpha),
             xlab = 'Homologue Specific Relative CN')
    abline(v = (alpha*c(0:6) + 1-alpha)/mean(c(a1, a2)), col = 'green')
    #  delta2 = mean(diff(h$mids[h$counts >0]))
    #  bs = min(h$mids[h$counts>0]) - delta2*c(0:2)
    #  (1-alphas)/alphas*delta2
  }
  
  #Assume first solution is correct (what happens to the different cases)
  par(mfcol = c(2,3))
  alphas = c(.28, .39, .64)
  for (i in 1:3){
    alpha = alphas[i]
    a1 = c(0, 1, 1, 2, 2) + i - 1
    a2 = c(2, 1, 2, 2, 3) + i - 1
    a1 = alpha*a1 +(1-alpha)*1
    a2 = alpha*a2 +(1-alpha)*1
    h = hist(c(a1, a2), breaks = 100, xlim = c(0, 3.5))
    abline(v = alpha*c(0:6) + 1-alpha, col = 'green') 
    alphaS = mean(diff(h$mids[h$counts >0])) #Change
    S = alphaS/alphas
    b0 = min(h$mids[h$counts>0]) - alphaS*c(0:2) #Change
    oneMinusAlpha = b0/S
    alpha = 1-oneMinusAlpha
    frame()
    legend("topleft", legend = c(paste('alphaS =', alphaS),
                                 paste('S =', S),
                                 paste('b0 =', b0),
                                 paste('oneMinusAlpha =',oneMinusAlpha),
                                 paste('alpha =', alpha)),
           bty = 'n')
  }
  
  alpha = c(.28, .4, .62)
  tau = c(2.94, 4.91, 6.87)
  ybar = tau/2
  delta2 = mean(diff(h$mids[h$counts>0]))
  b = min(h$mids[h$counts>0])
  b = b -delta2*c(0, 1, 2)
  delta2/alpha
  b/(1-alpha)
  
  
  
  
  par(mfcol = c(2,3))
  alphas = c(.28, .39, .64)
  for (i in 1:3){
    f = .7
    alpha = alphas[i]
    a1 = c(0, 1, 1, 2, 2) + i - 1
    a2 = c(2, 1, 2, 2, 3) + i - 1
    a1 = f*(alpha*a1 +(1-alpha)*1)
    a2 = f*(alpha*a2 +(1-alpha)*1)
    h = hist(c(a1, a2), breaks = 100, xlim = c(0, 3.5), main  = paste('alpha =', alpha),
             xlab = 'Homologue Specific CN')
    abline(v = alpha*c(0:6) + 1-alpha, col = 'green') 
    r = c(a1, a2)/mean(c(a1, a2))
    h = hist(r, breaks = 100, xlim = c(0, 1.5), main  = paste('alpha =', alpha),
             xlab = 'Homologue Specific Relative CN')
    abline(v = (alpha*c(0:6) + 1-alpha)/mean(c(a1, a2)), col = 'green')
    #  delta2 = mean(diff(h$mids[h$counts >0]))
    #  bs = min(h$mids[h$counts>0]) - delta2*c(0:2)
    #  (1-alphas)/alphas*delta2
  }
  
  a1 = f*(alpha*A1 +(1-alpha)*1)
  a2 = f*(alpha*A2 +(1-alpha)*1)
  
  
  
  
  ################################################################
  ################################################################
  #
  # Birdseed with Affymetrix Power Tools
  #
  ################################################################
  ################################################################
  
  
  ################################################################
  ##  Run APT from anywhere
  
  
  #Run apt-probeset-summarize from RESPONSIFY/HAPSEG
  apt-probeset-summarize --cdf-file  ~/RESPONSIFY/AROMA/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --analysis rma-bg,quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir ~/RESPONSIFY/ABSOLUTE/APT/ ~/RESPONSIFY/AROMA/rawData/responsify/GenomeWideSNP_6/*.CEL
  
  #This file (or similar) will be your snp.fn argument to HAPSEG:
  rma-bg.quant-norm.pm-only.med-polish.expr.summary.txt 
  
  #Next you will need to run birdseed:
  apt-probeset-genotype -o BirdseedOutput -c ~/RESPONSIFY/AROMA/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --set-gender-method cn-probe-chrXY-ratio --chrX-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrXprobes --chrY-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrYprobes --special-snps ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.specialSNPs --read-models-birdseed ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.v2.6.birdseed.models -a birdseed-dev --write-models ~/RESPONSIFY/AROMA/rawData/responsify/GenomeWideSNP_6/*.CEL
  
  #These will be input to HAPSEG:
  birdseed-dev.calls.txt will be calls.fn
  birdseed-dev.snp-models.txt which will be clusters.fn  
  
  #Also let
  snp.file.parser=AptSnpFileParser
  cluster.file.parser=BirdseedClustersFileParser
  
  ################################################################
  ##  Run APT from anywhere ATTEMPT 2
  
  
  #Run apt-probeset-summarize from RESPONSIFY/HAPSEG
  apt-probeset-summarize --cdf-file  ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir ~/RESPONSIFY/ABSOLUTE/APT2/ ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL
  
  #This file (or similar) will be your snp.fn argument to HAPSEG:
  quant-norm.pm-only.med-polish.expr.summary.txt 
  
  #Next you will need to run birdseed:
  apt-probeset-genotype -o BirdseedOutput -c ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --set-gender-method cn-probe-chrXY-ratio --chrX-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrXprobes --chrY-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrYprobes --special-snps ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.specialSNPs --read-models-birdseed ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.v2.6.birdseed.models -a birdseed-v1 --write-models ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL
  
  #These will be input to HAPSEG:
  birdseed-v1.calls.txt will be calls.fn
  birdseed-v1.snp-models.txt which will be clusters.fn  
  
  #Also let
  snp.file.parser=AptSnpFileParser
  cluster.file.parser=BirdseedClustersFileParser
  
  
  ################################################################
  ################################################################
  #
  # HAPSEG from Birdseed
  #
  ################################################################
  ################################################################
  prefix.out2 = paste(getwd(), '/ABSOLUTE/HAPSEGbirdseed/', sep = '')
  library(HAPSEG)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
  info$sample.name = paste(info$ind, info$sample, sep = '.')
  
  DoHapseg <- function(i, info) {
    myarray = info$arrname[i]
    mysample = info$sample.name[i]
    normal = info$tissue[i] == 'N'
    registerDoSEQ()
    plate.name <- mysample
    clusters.fn <- paste(getwd(), "/ABSOLUTE/APT/BirdseedOutput/birdseed-dev.snp-models.txt", sep = '')  
    phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
    force.diploid = normal
    d = file.path(getwd(), 'ABSOLUTE', 'APT')
    snp.fn <- file.path(d,"rma-bg.quant-norm.pm-only.med-polish.expr.summary.txt")
    calls.fn = file.path(d, 'BirdeseedOutput', 'birdseed-dev.calls.txt')
    mn.sample = NULL
    if (!normal) {
      mn.sample = info$arrname[info$ind == info$ind[i]+100]
      if (length(mn.sample) == 0) mn.sample = NULL
    }
    results.dir <- paste(prefix.out2, myarray, sep = '')
    if (!file.exists(results.dir)) {
      dir.create(results.dir, recursive=TRUE)
    }
    print(paste("Starting array", mysample, "at", results.dir))
    log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEGbirdseed', 'log')
    if (!file.exists(log.dir)) {
      dir.create(log.dir, recursive=TRUE)
    }
    sink(file=file.path(log.dir, paste(myarray, ".hapseg.out.txt", sep="")))
    myRunHapSeg(plate.name, myarray, seg.fn=NULL, snp.fn, genome = 'hg19', results.dir, platform = "SNP_6.0", 
                use.pop='CEPH', 
                impute.gt=TRUE, plot.segfit=TRUE, merge.small=TRUE, merge.close=TRUE, min.seg.size=5, 
                normal=normal, 
                out.p=0.001, seg.merge.thresh=1e-10, use.normal = TRUE, adj.atten = FALSE, 
                phased.bgl.dir, calibrate.data=TRUE, 
                drop.x=FALSE, verbose=TRUE, snp.file.parser=AptSnpFileParser,
                calls.fn = calls.fn, force.diploid = force.diploid, clusters.fn = clusters.fn,
                clusters.file.parser=BirdseedClustersFileParser)
    sink()
  }
  
  foreach (i=1:86, .combine=c) %dopar% {
    DoHapseg(i=i, info = info)
  }  
  
  