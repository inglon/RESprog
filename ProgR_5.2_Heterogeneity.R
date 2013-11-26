################################################################
################################################################
# Script: ProgR_5.1_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  23/05/2013
# R version: R 3.0
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

prefix.raw = paste(getwd(), "/AROMA/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/AROMA/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/ABSOLUTE/HAPSEG/", sep='')


date = format(Sys.Date())




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

######Find segments with 2m1 and 2m0

pdf(paste(prefix.out, 'Find_segments_with_2m1_2m0_', date, '.pdf', sep=''), paper = 'a4r',
    width = 11, height = 8) 
par(mfrow=c(2,2), mar = c(2,2,3,0)+.1)
for (i in 1:56){
ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                          paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))
c2 = subset(ta, Cn == 2 & !is.na(imba))
plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
points(imba~log2, data = c2, col = 'red')
tmp = sort(kmeans(c2$imba, 2)$centers)
het = tmp[1]
loh = tmp[2]
mid = (het + loh)/2  

loh.segs = subset(c2, imba>mid)
if (nrow(loh.segs)>1){
den = density(loh.segs$imba)
loh.r = range(den$x[which(den$y>max(den$y)*2/3)])
loh.segs = subset(c2, imba>loh.r[1] & imba<loh.r[2], select = c('Chromosome','Start','End'))
points(imba~log2, data=c2, subset = imba>loh.r[1] & imba<loh.r[2], col = 'green')
} else {
  loh.segs = subset(c2, select = c('Chromosome','Start','End'))
  points(imba~log2, data=c2, subset = imba>mid, col = 'green')
}

het.segs = subset(c2, imba<mid)
den = density(het.segs$imba)
het.r = range(den$x[which(den$y>max(den$y)/2)])
het.segs = subset(c2, imba>het.r[1] & imba<het.r[2], select = c('Chromosome','Start','End'))
points(imba~log2, data=c2, subset = imba>het.r[1] & imba<het.r[2], col = 'blue')

d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
dir.create(d, recursive= T, showWarnings = F)
write.table(het.segs, file.path(d, 'Hetsegs.txt'), row.names = F, quote = F, sep = '\t')
write.table(loh.segs, file.path(d, 'Homsegs.txt'), row.names = F, quote = F, sep = '\t')
}
dev.off()
     
######Find allele intensities

for (i in 1:56){   
library(aroma.affymetrix)
dataSet <- "responsify";
tags <- "ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2";
chipType <- "GenomeWideSNP_6";
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
df <- getFile(ds, which(getNames(ds) == info$arrname2[i]));
bs <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
bf <- getFile(bs, which(getNames(ds) == info$arrname2[i]));
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


######Find calibration constants
i = 24
i = 21
i = 26

library(IRanges)
pdf(paste(prefix.out, 'Find_calibration_constants_', date, '.pdf', sep=''), paper = 'a4') 
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  homs = read.delim(file.path(d, 'Homsegs.txt'), stringsAsFactors = F)
  hets = read.delim(file.path(d, 'Hetsegs.txt'), stringsAsFactors = F)
  data = read.delim(file.path(d, 'Total,fracB.txt'), stringsAsFactors = F)
  
  nc = nchar(homs$Chromosome, type = 'c')
  homs$Chromosome = substr(homs$Chromosome, 4, nc)
  homs$Chromosome[homs$Chromosome == 'X'] = '23'
  homs$Chromosome[homs$Chromosome == 'Y'] = '24'
  homs$Chromosome[homs$Chromosome == 'M'] = '25'
  homs$Chromosome = as.integer(homs$Chromosome)
  
  nc = nchar(hets$Chromosome, type = 'c')
  hets$Chromosome = substr(hets$Chromosome, 4, nc)
  hets$Chromosome[hets$Chromosome == 'X'] = '23'
  hets$Chromosome[hets$Chromosome == 'Y'] = '24'
  hets$Chromosome[hets$Chromosome == 'M'] = '25'
  hets$Chromosome = as.integer(hets$Chromosome)
  
  subd = subset(data, !is.na(chromosome) & !is.na(position) & !is.na(fracB) & !is.na(total) & !is.na(unitNames))
  subd$cn2 = 'no'
  subd$chromosome = factor(subd$chromosome)
  homs$Chromosome = factor(homs$Chromosome, levels = levels(subd$chromosome))
  h = RangedData(IRanges(start = homs$Start, end = homs$End),  
                 space = homs$Chromosome, universe = 'hg19') 
  d = RangedData(IRanges(start = subd$position, width = 1),  
                 space = subd$chromosome, universe = 'hg19')
  hits = findOverlaps(query = h, subject = d)
  for (ch in levels(subd$chromosome)){
    myhits = subjectHits(hits[[ch]])
    ix = which(subd$chromosome == ch)
    subd$cn2[ix[myhits]] = 'hom' 
  }
  
  hets$Chromosome = factor(hets$Chromosome, levels = levels(subd$chromosome))
  h = RangedData(IRanges(start = hets$Start, end = hets$End),  
                 space = hets$Chromosome, universe = 'hg19') 
  hits = findOverlaps(query = h, subject = d)
  for (ch in levels(subd$chromosome)){
    myhits = subjectHits(hits[[ch]])
    ix = which(subd$chromosome == ch)
    subd$cn2[ix[myhits]] = 'het' 
  }
  
  subd = subset(subd, cn2 %in% c('hom','het'))
  subd$baf = pmax(subd$fracB, 0)
  subd$baf = pmin(subd$fracB, 1)
  
  subd$a = subd$total*(1-subd$baf)
  subd$b = subd$total*subd$baf
  
  par(mfrow=c(2,2))
  #A0
  x = subd$a[subd$baf>.5 & subd$cn2 == 'hom']
  q = quantile(x, c(.05, .9))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'A0: A intensities at BB')
    A0 = den$x[which.max(den$y)]
  } else {
    A0 = ifelse(length(x)>0, x, 0)
  }                    
  
  #B0
  x = subd$b[subd$baf<.5 & subd$cn2 == 'hom']
  q = quantile(x, c(.05, .9))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'B0: B intensities at AA')
    B0 = den$x[which.max(den$y)]
  } else {
    B0 = ifelse(length(x)>0, x, 0)
  }
  
  #A1
  x = subd$a[subd$baf>.3 & subd$baf<0.7 & subd$cn2 == 'het']
  q = quantile(x, c(.05, .95))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'A1: A intensities at AB')
    A1 = den$x[which.max(den$y)]
  } else {
    A1 = ifelse(length(x)>0, x, 1)
    
  }
  
  #B1
  x = subd$b[subd$baf>.3 & subd$baf<0.7 & subd$cn2 == 'het']
  q = quantile(x, c(.05, .95))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'B1: B intensities at AB')
    B1 = den$x[which.max(den$y)]
  } else {
    B1 = ifelse(length(x)>0, x, 1)
    
  }
  
  title(main = info$sample.name[i], outer = T, line = -1)
  
  ###Standardize the A and B intensities and write to file for HAPSEG
  
  data$baf = pmax(data$fracB, 0)
  data$baf = pmin(data$fracB, 1)
  data$a = data$total*(1-data$baf)
  data$b = data$total*(data$baf)
  data$A = (data$a-A0)/(A1-A0)
  data$B = (data$b-B0)/(B1-B0)
  
  dat = subset(data, select = c('A','B'))
  rownames(dat) = data$unitNames
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  #write.table(data, file.path(d, 'snpfile.txt'), quote = F, sep = '\t')
  save(dat, file = file.path(d, 'snpfile.RData'))
}
dev.off()





###### Segment Logratios for calibrated data

library(aroma.affymetrix)
tot = NULL
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  load(file.path(d, 'snpfile.RData'))
  tot = cbind(tot, apply(dat, 1, sum))
  rm(dat)
}
den = apply(tot, 1, median, na.rm = T)

dataSet <- "responsify";
tags <- "ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2";
chipType <- "GenomeWideSNP_6";
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
df <- getFile(ds, which(getNames(ds) == info$arrname2[1]));
ugp <- getAromaUgpFile(df)
chromosome <- ugp[, 1, drop = TRUE]
position <- ugp[, 2, drop = TRUE]

data = data.frame(chromosome=chromosome, position = position, 
                  stringsAsFactors = F)

library(PSCBS)
library(IRanges)
for (i in 1:56){
  print(i)
  #Prepare data for this sample
  fitname = paste('fit', i, sep = '')
  filename <- paste(getwd(), '/output/SNP/CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  seg = getSegments(fit)
  seg$Sample = info$myarray[i]
  seg = data.frame(Sample = seg$Sample, Chromosome = seg$chromosome, "Start Position" = seg$tcnStart,
                   "End Position" = seg$tcnEnd, "Num markers" = seg$tcnNbrOfSNPs, 
                   check.names = F)
  seg = subset(seg, !is.na(Chromosome) & Chromosome != 25)
  seg$Chromosome = factor(seg$Chromosome, levels = as.character(1:24))
  seg$cn = NA
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  load(file.path(d, 'snpfile.RData'))
  data$cn = (apply(dat, 1, sum)/den)
  ix = !is.na(data$chromosome) & !is.na(data$position) & !is.na(data$cn) & data$chromosome != 25
  h = RangedData(IRanges(start = seg[,'Start Position'], end = seg[,'End Position']),  
                 space = seg$Chromosome, universe = 'hg19') 
  d = RangedData(IRanges(start = data$position[ix], width = 1),  
                 space = data$chromosome[ix], universe = 'hg19', cn = data$cn[ix])
  hits = findOverlaps(query = h, subject = d)
  
  for (ch in levels(seg$Chromosome)){
    cn = d$cn[space(d) == as.numeric(ch)]
    seg$cn[seg$Chromosome == ch & seg[,'Num markers']>0] = tapply(cn[subjectHits(hits[[ch]])], queryHits(hits[[ch]]), 
                                          median, na.rm = T)
  }
  seg$cn = log2(seg$cn)
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i], 'Segdat.txt')
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
#source(paste(getwd(), '/programs/Hapseg.r', sep = ''))
library(foreach)
library(doMC)
registerDoMC(20)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(info$ind, info$sample, sep = '.')

myRunHapSeg = function (plate.name, array.name, seg.fn, snp.fn, genome.build, 
                        results.dir, platform, use.pop, impute.gt, plot.segfit, merge.small, 
                        merge.close, min.seg.size, normal, out.p, seg.merge.thresh, 
                        use.normal, adj.atten, phased.bgl.dir, force.diploid = normal, 
                        drop.x = FALSE, drop.y = TRUE, calls.fn = NULL, mn.sample = NULL, 
                        out.file = NULL, calibrate.data = NULL, clusters.fn = NULL, 
                        prev.theta.fn = NULL, snp.file.parser = DefaultSnpFileParser, 
                        clusters.file.parser = DefaultClustersFileParser, verbose = FALSE) 
{
  if (!HAPSEG:::CheckGenomeBuild(genome.build)) {
    stop("Unsupported genome build: ", genome.build)
  }
  if ((!impute.gt) && (verbose)) {
    print("Note: impute.gt is set to FALSE. It is recommended that this be set to TRUE")
  }
  if (!identical(adj.atten, FALSE)) {
    print("Note: The adj.atten argument has been disabled and will be removed in a future release")
  }
  tmp.dir = HAPSEG:::CreateTmpDir(results.dir)
  on.exit(try(unlink(tmp.dir, recursive = TRUE), silent = TRUE), 
          add = TRUE)
  iams.res = HAPSEG:::InitAndMergeSmall(array.name, genome.build, use.pop, 
                                        use.normal, normal, impute.gt, platform, seg.fn, snp.fn, 
                                        drop.x, drop.y, calls.fn, mn.sample, min.seg.size, merge.close, 
                                        out.p, merge.small, force.diploid, calibrate.data, clusters.fn, 
                                        prev.theta.fn, snp.file.parser = snp.file.parser, clusters.file.parser = clusters.file.parser, 
                                        verbose = verbose)
  if (merge.close == TRUE) {
    mcaf.res = HAPSEG:::MergeCloseAndFit(iams.res, out.p, seg.merge.thresh, 
                                         impute.gt, platform, force.diploid = force.diploid, 
                                         verbose = verbose)
  }
  else {
    if (verbose) {
      print("Merged Loci: nulled ")
    }
    seg.dat[["merged.loci"]] = NA
    seg.dat[["final.merge.prob"]] = matrix(NA, ncol = 2, 
                                           nrow = length(iams.res[["h.d"]]) - 1)
    seg.dat[["seg.expected.phase"]] = iams.res[["em.fit"]][["seg.expected.phase"]]
    mcaf.res = iams.res
    mcaf.res[["seg.dat"]] = seg.dat
  }
  if (impute.gt == TRUE) {
    gt.res = HAPSEG:::ImputeGt(mcaf.res, platform, tmp.dir, plate.name, 
                               out.p, phased.bgl.dir, verbose = verbose)
    gt.res = HAPSEG:::PostImputeSegFit(gt.res, out.p, platform, force.diploid = force.diploid, 
                                       verbose = verbose)
  }
  else {
    gt.res = mcaf.res
  }
  seg.dat = HAPSEG:::FinishHapsegAndSave(gt.res, plate.name, normal, 
                                         array.name, results.dir, platform, out.file)
  if (plot.segfit) {
    if (nrow(seg.dat$merged.loci) == 0) seg.dat$merged.loci = NA
    HAPSEG:::DoPlots(seg.dat, results.dir, platform, verbose = verbose)
  }
  return(TRUE)
}

DoHapseg <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- paste(prefix.out, mysample, sep = '')
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- paste(prefix.out, 'log', sep = '')
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
  DoHapseg(myarray = info$myarray[i], mysample = info$sample.name[i])
}


#Discard:
DoHapseg2 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  results.dir <- paste(prefix.out2, mysample, sep = '')
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- paste(prefix.out, 'log', sep = '')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn=NULL, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=F, merge.close=T, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-3, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}



foreach (i=c(21, 15, 24, 26), .combine=c) %dopar% {
  DoHapseg2(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}


#Discard:
for (arr in 28:56){
  sample.name <- info$sample.name[arr]
  mydir = paste(prefix.out, sample.name, sep ='')
  for (j in 1:22){
    subdir = file.path(mydir,paste('chr', j, sep = ''))
    allfiles = list.files(subdir)
    tmp = unlist(lapply(strsplit(allfiles, split = '_'), function(x){x[3]}))
    nc = nchar(tmp, type = 'c')
    tmp = substr(tmp, 1, nc-4)
    ix = sort(as.numeric(tmp))
    conseq = (min(ix)):(max(ix))
    gap = setdiff(conseq, ix)
    if (length(gap)>0){
      delix = min(ix):(min(gap)-1)
      delfiles = paste('HAPSEG_SEG_', delix, '.jpg', sep = '')
      file.remove(file.path(subdir, delfiles))          
    }
  }
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
# Absolute TOTAL copy numbers in segments: 
# This probably did not find the mutation file, so not useful.
# 
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
prefix.ABS = paste(getwd(), '/ABSOLUTE/RunAbsoluteTCN/', sep = '')
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname


DoAbsoluteTCN <- function(i) {
  registerDoSEQ()
  sample.name = paste(i, info$sample[i], sep = '.')
  fitname = paste('fit', i, sep = '')
  seg.dat.fn = paste(prefix.ABS, 'seg.dat.files/', fitname, 'tcn.txt', sep = '');
  results.dir = prefix.ABS
  maf.fn = paste(getwd(), 
                 '/rawData/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_',
                 sample.name, '.maf', sep = '')
  print(paste("Starting myarray", sample.name, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsoluteTCN/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  
  sink(file=file.path(log.dir, paste(i, '.', sample.name, ".abs.out.txt", sep="")))
  RunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 2, max.ploidy = 4, 
              primary.disease = 'Breast Cancer', platform = "SNP_6.0", sample.name = sample.name, 
              results.dir = results.dir,
              max.as.seg.count = 1500, max.non.clonal = 0.5, max.neg.genome = 0, copy_num_type = 'total', 
              maf.fn=maf.fn, min.mut.af=NULL, output.fn.base=sample.name,
              verbose=T)
  #Note, min and max ploidy might have to be changed between samples
  #As would max.non.clonal and the rest...
  sink()
}

foreach (i = 1:56, .combine=c) %dopar% {
  DoAbsoluteTCN(i)
}

###Follow up sample 24
sample.name = paste(i, info$sample[i], sep = '.')
CreateReviewObject(obj.name=sample.name, absolute.files=paste(prefix.ABS,'24.4742.ABSOLUTE.RData', sep = ''), 
                   prefix.ABS, "total", verbose=TRUE)
calls.path = paste(prefix.ABS, sample.name,".PP-calls_tab.txt", sep = '')
modes.path = paste(prefix.ABS, sample.name, ".PP-modes.data.RData", sep = '')
output.path = paste(prefix.ABS, sample.name, "abs_extract", sep = '')
ExtractReviewedResults(calls.path, "IL", modes.path, output.path, "absolute", "total")

#Just manually checking what the output looks like
absdat = read.delim(paste(output.path,'/reviewed/SEG_MAF/24.4742.segtab.txt', sep = ''))
hist(absdat$cancer_cell_frac[!(absdat$cancer_cell_frac %in% c(0,1))], breaks = 100) #Cancer cell fractions
na.omit(absdat$cancer_cell_frac[!(absdat$cancer_cell_frac %in% c(0,1))])
[1] 0.89 0.59 0.66 0.64 0.68 0.49 0.49 0.56 0.63 0.52 0.57
na.omit(absdat$subclonal[!(absdat$subclonal %in% c(0,1))])




################################################################
################################################################
#
# Run ABSOLUTE
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
  seg.dat.fn <- paste(prefix.out, mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  RunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


foreach (i = c(2, 3, 12:14,16, 20, 22:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}


foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}


d = paste(getwd(), '/ABSOLUTE/RunAbsolute', sep = '')
obj.name = 'ALL'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                    info$sample.name[c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56)], 
                                                           '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)



#Not done yet:
obj.name = 'ALL'
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "test", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "default", copy_num_type="allelic")



################################################################
################################################################
#
# Check absolute solution
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
################################################################
######Find clusters of segments which we like to check

for (i in 1:56){
  
  i = 56
  ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                            paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  
  d0.imba = diff(range(ta$imba, na.rm = T))/200
  d0.log = diff(range(ta$log2, na.rm = T))/1000
  out = NULL
  
  n = 4 ###MANUALLY INSERT THE NUMBER OF CLUSTERS YOU LIKE HERE:
  pos = locator(n)
  
      ###FOR EACH CLUSTER (CHANGE c BETWEEN 1, 2, 3, ..., n)
  c = 4
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
######Find HSCRs for the chosen clusters
library(IRanges)
col.transp <- function(col1, alpha=0.5){
  tmp=as.list(t(col2rgb(col1)))
  tmp$alpha=alpha*255
  tmp$maxColorValue=255
  do.call('rgb',tmp)
}
cols = c(2:7, 'gold' , colors()[c(96)])
for (c in 7:length(cols)) cols[c] = col.transp(cols[c])


obj.name = 'ALL'
d = paste(getwd(), '/ABSOLUTE/RunAbsolute', sep = '')
modes.fn = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))

ABSOLUTE:::set_allelic_funcs()
load(modes.fn)

for (i in c(15, 21)){
  segobj = segobj.list[[info$sample.name[i]]]
  if (!is.null(segobj$maf.fn)) if (file.exists(segobj$maf.fn)){
    
    ### Load data for this sample
    ###
    
    d = file.path(getwd(), 'ABSOLUTE', 'Choose_solution')
    pdf(paste(d, '/', info$sample.name[i],'_Solutions_output', '.pdf', sep=''), 
        paper = 'a4', height = 12, width = 8)   
    d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  subsegs = read.delim(file.path(d,'Subsegs.txt'))
  homsegs = read.delim(file.path(d, 'Homsegs.txt'))
  homsegs$cluster = 'hom'
    hetsegs = read.delim(file.path(d, 'Hetsegs.txt'))
    hetsegs$cluster = 'het'
    cl = rbind(subsegs, homsegs, hetsegs)
    cl = subset(cl, !(Chromosome %in% c('chrX', 'chrY')))
    cl$Chromosome = as.character(cl$Chromosome)
    nc = nchar(cl$Chromosome, type = 'c')
    cl$Chromosome = as.numeric(substr(cl$Chromosome, 4, nc))
    cl$Chromosome = factor(cl$Chromosome, levels = as.character(1:22))
    
    
    
    nsol = nrow(segobj$mode.res$mode.tab)
    for (j in 1:nsol){
      
      ##Extact absolute solution data
      segtmp = segobj
      segtmp$mode.res = ABSOLUTE:::ReorderModeRes(segtmp$mode.res, j)
      purity = round(segtmp$mode.res$mode.tab[,'alpha'], 2)
      ploidy = round(segtmp$mode.res$mode.tab[,'tau'], 2)  
      frac.het = segtmp$mode.res$mode.tab[,'frac.het']  
      at = segtmp$mode.res$mode.tab[1, "AT"]
      ab = as.data.frame(GetAbsSegDat(segtmp, at))
      ab$Chromosome = factor(ab$Chromosome, levels = as.character(1:22))
      
      ##Index cluster overlaps with absolute segments
      for (c in seq_along(unique(cl$cluster))){
        cc = subset(cl, cluster == unique(cl$cluster)[c])
        ab[,paste('cluster', unique(cl$cluster)[c])] = 0
        for (ch in seq_along(levels(cc$Chromosome))){
          ab.ix = which(ab$Chromosome == levels(cc$Chromosome)[ch])
          cc.ix = which(cc$Chromosome == levels(cc$Chromosome)[ch])
          IRab = IRanges(start = ab$Start.bp[ab.ix], end = ab$End.bp[ab.ix])
          IRcc = IRanges(start = cc$Start[cc.ix], end = cc$End[cc.ix])
          if (length(IRcc)>0 & length(IRab)>0){
            cov = coverage(c(IRab, IRcc)) - 1
            ab[ab.ix,paste('cluster', unique(cl$cluster)[c])] = viewSums(Views(cov, IRab))
          } 
        }  
      }
      
      
      ### Plots of HAPSEG output
      ###
      if (j == 1) {
        par(mfrow = c(3,2))
        #imba by logratio from hscrs
        imba = 2*abs(ab$hscr.a1/(ab$hscr.a1+ab$hscr.a2) - 0.5) 
        lr = log(ab$copy.ratio)
        size = ab$W*100
        Wmod = lm(W~length, data = ab)
        ran = quantile(lr, prob = c(.05, .95))
        ran = ran + c(-1, 1)*diff(ran)/10
        plot(lr, imba, xlim = ran, cex = size, xlab = 'log(Total Copy Number)',
             ylab = 'Alleleic imbalance 2*|BAF - 0.5|', main = 'HAPSEG Allelic Imbalance')
        ix = which(ab[,'cluster het']>0)
        points(lr[ix], imba[ix], col =cols[8], pch = 16, cex = size[ix])
        try({ix = which(ab[,'cluster hom']>0)
        points(lr[ix], imba[ix], col =cols[7], pch = 16, cex = size[ix])}
            , silent = T)
        for (c in seq_along(unique(cl$cluster))){
          clust = unique(cl$cluster)[c]
          if (!(clust %in% c('het','hom'))){
            ix = which(ab[,paste('cluster', clust)]>0)
            clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))*100
            points(lr[ix], imba[ix], col =cols[c], 
                   pch = 16, cex = clust.W[ix])
          }
        }
        clusters = unique(cl$cluster)
        nc = length(clusters)
        clusters = clusters[c(nc-1, nc, 1:(nc-2))]
        lcols = cols[c(7, 8, 1:(nc-2))]
        legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                paste('Cluster', 1:(nc-2))), bty = 'n')
        
        
        
        #hscr.a1 vs hscr.a2
        y = log(ab$hscr.a1)
        x = log(ab$hscr.a2)
        ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
        size = ab$W*100
        ranx = quantile(x[ix], prob = c(.05, .95))
        ranx = ranx + c(-1, 1)*diff(ranx)/10
        rany = quantile(y[ix], prob = c(.05, .95))
        rany = rany + c(-1, 1)*diff(rany)/10
        plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Allele 2 logratio',
             ylab = 'Allele 1 logratio', main = 'HAPSEG Copy Ratios')
        ix = which(ab[,'cluster het']>0)
        points(x[ix], y[ix], col =cols[8], pch = 16, cex = size[ix])
        try({ix = which(ab[,'cluster hom']>0)
        points(x[ix], y[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
        for (c in seq_along(unique(cl$cluster))){
          clust = unique(cl$cluster)[c]
          if (!(clust %in% c('het','hom'))){
            ix = which(ab[,paste('cluster', clust)]>0)
            clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))*100
            points(x[ix], y[ix], col =cols[c], pch = 16,, cex = clust.W[ix])
          }
        }
        clusters = unique(cl$cluster)
        nc = length(clusters)
        clusters = clusters[c(nc-1, nc, 1:(nc-2))]
        lcols = cols[c(7, 8, 1:(nc-2))]
        legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                                paste('Cluster', 1:(nc-2))), bty = 'n')
        
        
        ### Plots of histograms
        myab = ab
        myab$color.by = 8
        Wmod = lm(W~length, data = ab)
        nclust = length(unique(cl$cluster))
        colvec = c(as.numeric(cols[c(1:(nclust-2))]), 6, 7)
        for (c in seq_along(unique(cl$cluster))){
          clust = unique(cl$cluster)[c]
          mycl = myab[,paste('cluster', clust)]
          if (clust == 'het') clust = nclust - 1
          if (clust == 'hom') clust = nclust
          
          clust.W = predict(Wmod, data.frame(length=mycl))
          clust.W[mycl == 0] = 0
          myab$W = pmax(0, myab$W - clust.W)
          ix = mycl>0
          toadd = myab[ix,]
          toadd$W = clust.W[ix]
          toadd$color.by= colvec[c]
          if (clust == nclust-1) toadd$color.by = colvec[nclust-1]
          if (clust == nclust) toadd$color.by = colvec[nclust]
          myab = rbind(myab, toadd)
        }
        x = c(myab$hscr.a1, myab$hscr.a2)
        l = rep(myab$W, 2)
        colseg = rep(myab$color.by, 2)
        
        ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, 
                ylab = '', 
                data.whiskers = F,
                xlab = 'Pooled haplytype copy values', xlim = c(0, 3))
        grid()
        clusters = unique(cl$cluster)
        nc = length(clusters)
        clusters = clusters[c(nc-1, nc, 1:(nc-2))]
        lcols = colvec[c(nc, nc-1, 1:(nc-2))]
        legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                             paste('Cluster', 1:(nc-2))), bty = 'n')
        
        x = myab$hscr.a1
        l = myab$W
        colseg = myab$color.by
        
        ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                                  data.whiskers = F,
                                  xlab = 'Allele 1 copy values', xlim = c(0, 3))
        grid()
        
        frame()
        x = myab$hscr.a2
        l = myab$W
        colseg = myab$color.by
        
        ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                                  data.whiskers = F,
                                  xlab = 'Allele 2 copy values', xlim = c(0, 3))
        grid() 
        
      }
      
      
      ### Plot diagnostics for each ABSOLUTE solution
      ###
      
       
      
      #imba by logratio from expected
      imba = 2*abs(ab$expected.a1/(ab$expected.a1+ab$expected.a2) - 0.5) 
      lr = log((ab$expected.a1+ab$expected.a2)/2)
      size = ab$W*100
      Wmod = lm(W~length, data = ab)
      ran = quantile(lr, prob = c(.05, .95))
      ran = ran + c(-1, 1)*diff(ran)/10
      plot(lr, imba, xlim = ran, cex = size[ix], xlab = 'log(ABSOLUTE CN)',
           ylab = 'ABSOLUTE Allelic Imbalance 2*|BAF - 0.5|')
      ix = which(ab[,'cluster het']>0)
      points(lr[ix], imba[ix], col =cols[8], pch = 16, cex = size[ix])
      try({ix = which(ab[,'cluster hom']>0)
      points(lr[ix], imba[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
      for (c in seq_along(unique(cl$cluster))){
        clust = unique(cl$cluster)[c]
        if (!(clust %in% c('het','hom'))){
          clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))*100
          ix = which(ab[,paste('cluster', clust)]>0)
          points(lr[ix], imba[ix], col =cols[c], pch = 16, 
                 cex = clust.W[ix])
        }
      }
      clusters = unique(cl$cluster)
      nc = length(clusters)
      clusters = clusters[c(nc-1, nc, 1:(nc-2))]
      lcols = cols[c(7, 8, 1:(nc-2))]
      legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                              paste('Cluster', 1:(nc-2))), bty = 'n')
      
      title(main = paste('ABSOLUTE',j,'Allelic Imbalance'))
      
      
      
      #expected.a1 vs expected.a2
      y = ab$expected.a1
      x = ab$expected.a2
      ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
      size = ab$W*100
      ranx = quantile(x[ix], prob = c(.05, .95))
      ranx = ranx + c(-1, 1)*diff(ranx)/10
      rany = quantile(y[ix], prob = c(.05, .95))
      rany = rany + c(-1, 1)*diff(rany)/10
      rany[1] = min(0, rany[1])
      ranx[1] = min(0, ranx[1])
      rany[2] = max(2, rany[2])
      ranx[2] = max(2, ranx[2])
      plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Allele 2 ABSOLUTE CN',
           ylab = 'Allele 1 ABSOLUTE CN', type = 'n')
      at = ceiling(par('usr')[1L]):floor(par('usr')[2L])
      abline(v = at, col = "lightgray", lty = "dotted",
             lwd = par("lwd"))
      at = ceiling(rany[1]):floor(rany[2])
      abline(h = at, col = "lightgray", lty = "dotted",
             lwd = par("lwd"))
      ix = which(ab[,'cluster het']>0)
      points(x[ix], y[ix], col =cols[8], pch = 16, cex = size[ix])
      try({ix = which(ab[,'cluster hom']>0)
      points(x[ix], y[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
      for (c in seq_along(unique(cl$cluster))){
        clust = unique(cl$cluster)[c]
        if (!(clust %in% c('het','hom'))){
          ix = which(ab[,paste('cluster', clust)]>0)
          clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))*100
          points(x[ix], y[ix], col =cols[c], pch = 16, cex = clust.W[ix])
        }
      }
      points(0:2, 2:0, pch = 3, cex = 5, font = 3, col = cols[c(7,8,7)])
      clusters = unique(cl$cluster)
      nc = length(clusters)
      clusters = clusters[c(nc-1, nc, 1:(nc-2))]
      lcols = cols[c(7, 8, 1:(nc-2))]
      legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                              paste('Cluster', 1:(nc-2))), bty = 'n')
      mtext(paste('Ploidy ', ploidy, ', ', 'Purity ', purity, sep = ''), 
            line = 0.1, adj = 0)
      title(main = paste('ABSOLUTE',j,'Copy Numbers'))
      
       
     }
    dev.off()
  }
}




###Junk

for (i in 1:56){
  segobj = segobj.list[[info$sample.name[i]]]
  maf.fn = segobj$maf.fn
  print(maf.fn)
}
  #  ALL_130516_Filter_CBS_MinDepth_20_56.4868.maf
  d = file.path(getwd(), 'AROMA', 'rawData','responsify','Exome', 'MAF')
  if (!is.null(maf.fn)) segobj.list[[info$sample.name[i]]]$maf.fn = file.path(d,strsplit(maf.fn, split = '/')[[1]][12])
}


x = c(myab$expected.a1, myab$expected.a2)
l = rep(myab$W, 2)
colseg = rep(myab$color.by, 2)

ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Single copy values', xlim = c(0, 3))
grid()
clusters = unique(cl$cluster)
nc = length(clusters)
clusters = clusters[c(nc-1, nc, 1:(nc-2))]
lcols = colvec[c(nc, nc-1, 1:(nc-2))]
legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                     paste('Cluster', 1:(nc-2))), bty = 'n')




#xlim = range(c(ab$expected.a1, ab$expected.a2))
#par(mfrow = c(2, 3))
#for (c in seq_along(unique(cl$cluster))){
#  clust = unique(cl$cluster)[c]
#  ABSOLUTE:::PlotSeglenHist(c(ab$expected.a1, ab$expected.a2), rep(ab[,paste('cluster', clust)], 2),
#                            color.by = c(ab$subclonal.a1, ab$subclonal.a2), ylab = '', 
#                            xlab = 'Allelic copy number interpretation', xlim = xlim)
#  mtext(paste('Cluster', clust))
#}
#title(main = paste('Solution', j), outer = T, line = -1)
#title(main = paste('Ploidy ', ploidy, ', ', 'Purity ', purity, sep = ''), line = -1, 
#      outer = T, adj = 0)


  x = NULL
  w = 0
  mycols = NULL
  for (c in seq_along(unique(cl$cluster))){
    clust = unique(cl$cluster)[c]
    xtmp = c(ab$expected.a1, ab$expected.a2) 
    wtmp = rep(ab[,paste('cluster', clust)], 2)
    x = c(x, xtmp)
    w = c(w, wtmp)
    trycol = as.numeric(clust)+1
    if (is.na(trycol)) trycol = ifelse(clust == 'het', 1, 8)
    mycols = c(mycols, rep(trycol, length(xtmp)))
  }
  xlim = range(c(ab$expected.a1, ab$expected.a2))
  ABSOLUTE:::PlotSeglenHist(x, w,
                            color.by = cols, ylab = '', 
                            xlab = 'Allelic copy number interpretation', xlim = xlim)
  title(main = paste('Solution', nsol)) 

  
plot(ab$expected.a1, ab$hscr.a1, log = 'y') 
plot(ab$expected.a2, ab$hscr.a2, log = 'y')
  


#Mark subclonal
ix = which(ab$subclonal.a1 >0 | ab$subclonal.a2 > 0)
points(lr[ix], imba[ix], pch = 4, cex = size[ix])

#Mark 2m0
ix = which(ab[,'cluster hom']>0 & ((ab$modal.a1 == 0 & ab$modal.a2 == 2) | (ab$modal.a1 == 2 & ab$modal.a2 == 0)))
ab[ix, c('modal.a1','modal.a2')]
points(lr[ix], imba[ix], pch = 2, cex = size[ix])

#Mark 2m1
ix = which(ab[,'cluster hom']>0 & (ab$modal.a1 == 1 & ab$modal.a2 == 1))
ab[ix, c('modal.a1','modal.a2')]
points(lr[ix], imba[ix], pch = 5, cex = size[ix])




################################################################
######Print ABSOLUTE output for subset of solutions


d = paste(getwd(), '/ABSOLUTE/RunAbsolute', sep = '')
obj.name = '15_subset'
myCreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', info$sample.name[c(15)], 
                                                           '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE, myModes = c(4, 6, 11, 12, 13, 14))

obj_name = 'tmp'
myCreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', info$sample.name[c(15)], 
                                                             '.ABSOLUTE.RData', sep = ''), 
                     indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                     plot.modes = TRUE, verbose=TRUE, myModes = c(4, 6, 11, 12, 13, 14))



###############################################################

#ploidy check
i = 15
ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                          paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))

# This function estimates the median ploidy from a segmentation matrix of PICNIC results (integer copy number)
# by calculating the weighted median ploidy. The weights are equal to the segment length.
# The variabel seg needs to me a segmentation matrix with columns in the follwing order:
# Sample Name, Chromosome, Start Region, End Region, Number Probes, Copy Number Value.
# The matrix for seg must be a character matrix and mustn't be a dataframe.
# The variable x is should in include the unique sample names in the first column of seg .
# The function returns the median ploidy for each sample.

s = data.frame(Sample = i, Chromosome = ta$Chromosome, Start = ta$Start, End = ta$End,
               probes = ta$probes, cn = ta$Cn)
s = as.matrix(s)
pl = fun.ploidy(as.character(i), s) #Gave 3 for both 15 and 21


################################################################
################################################################
#
# Sample 15
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')

i = 15
d = paste(prefix.out, info$sample.name[i], sep = '')
list.files(d)
load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
mat = as.data.frame(seg.dat$allele.segs)
names(mat)[6:7] = c('A1','A2')
  
###Read selected segments
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
subsegs = read.delim(file.path(d,'Subsegs.txt'))
homsegs = read.delim(file.path(d, 'Homsegs.txt'))
homsegs$cluster = 'hom'
hetsegs = read.delim(file.path(d, 'Hetsegs.txt'))
hetsegs$cluster = 'het'
cl = rbind(subsegs, homsegs, hetsegs)
cl = subset(cl, !(Chromosome %in% c('chrX', 'chrY')))
cl$Chromosome = as.character(cl$Chromosome)
nc = nchar(cl$Chromosome, type = 'c')
cl$Chromosome = as.numeric(substr(cl$Chromosome, 4, nc))
cl$Chromosome = factor(cl$Chromosome, levels = as.character(1:22))

###Assign cluster overlaps and weights to segments
for (c in seq_along(unique(cl$cluster))){
  cc = subset(cl, cluster == unique(cl$cluster)[c])
  mat[,paste('cluster', unique(cl$cluster)[c])] = 0
  for (ch in seq_along(levels(cc$Chromosome))){
    mat.ix = which(mat$Chromosome == levels(cc$Chromosome)[ch])
    cc.ix = which(cc$Chromosome == levels(cc$Chromosome)[ch])
    IRmat = IRanges(start = mat$Start.bp[mat.ix], end = mat$End.bp[mat.ix])
    IRcc = IRanges(start = cc$Start[cc.ix], end = cc$End[cc.ix])
    if (length(IRcc)>0 & length(IRmat)>0){
      cov = coverage(c(IRmat, IRcc)) - 1
      mat[mat.ix,paste('cluster', unique(cl$cluster)[c])] = viewSums(Views(cov, IRmat))
    } 
  }  
}

#Transform
y = log(mat$A1)
x = log(mat$A2)
  
  



require(Rsolnp)
minfunc = function(pars, x, y, xclust, yclust){
  a = pars[1]
  b = pars[2]
  d = pars[3]
  xp = a*x + b*y
  yp = b*x + d*y
  sum(tapply(c(xp, yp), c(xclust, yclust), var))
}
eqfunc = function(pars, x, y, xclust, yclust){
  a = pars[1]
  b = pars[2]
  d = pars[3]
  a*d - b**2
}
pars = solnp(pars = c(1,1,1), fun = minfunc, 
             LB = c(0,0,0), x=x, y=y, 
             eqfun = eqfunc, eqB = 0)$pars
mat$A1new = alpha*mat$A1
mat$A2new = 1/alpha*mat$A2
sol = try(solnp(pars = c(mu0, r0), fun = minfunc, 
                LB = c(0,0), mymode = mymode, 
                o=fkcounts[as.numeric(names(fkcounts))>=atmin 
                           & as.numeric(names(fkcounts))<=atmax], 
                xmax = xmax), silent=T)

#Print new version
names(mat)[6:7] = c('A1.Seg.CN','A2.Seg.CN')
seg.dat$allele.segs = as.matrix(mat)
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled_segdat.RData', sep = '')))








y = (mat$A1)
x = (mat$A2)
ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
size = mat$length/sum(mat$length)*100
ranx = quantile(x[ix], prob = c(.05, .95))
ranx = ranx + c(-1, 1)*diff(ranx)/10
rany = quantile(y[ix], prob = c(.05, .95))
rany = rany + c(-1, 1)*diff(rany)/10
plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Haplotype 2 copy value',
     ylab = 'Haplotype 1 copy value', main = 'HAPSEG Copy Values')
x = c(mat$A1new, mat$A2new)
l = rep(size, 2)
ABSOLUTE:::PlotSeglenHist(x, l, color.by = NA, use.pal = NA, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Single copy values', xlim = c(0, 3))
grid()
x = c(mat$A1new)
l = (size)
ABSOLUTE:::PlotSeglenHist(x, l, color.by = NA, use.pal = NA, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Single copy values', xlim = c(0, 3))
grid()
x = c(mat$A2new)
l = (size)
ABSOLUTE:::PlotSeglenHist(x, l, color.by = NA, use.pal = NA, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Single copy values', xlim = c(0, 3))
grid()

  



